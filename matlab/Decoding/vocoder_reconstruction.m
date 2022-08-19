% Adapted from CNSP-Workshop 2021
% https://cnsp-workshop.github.io/website/index.html
% Author: Giovanni M. Di Liberto, Sophie Crowley
% Copyright 2021 - Giovanni Di Liberto
%                  Nathaniel Zuk
%                  Michael Crosse
%                  Aaron Nidiffer
%                  (see license file for details)
% Last update: 19 August 2022

clear all
close all

addpath libs/cnsp_utils
addpath libs/cnsp_utils/cnd
addpath libs/mTRF-Toolbox_v2/mtrf
addpath libs/NoiseTools
addpath eeglab2021.1


%% Set Parameters

% dataMainFolder = '../datasets/LalorNatSpeech/';
dataMainFolder = '../';
% dataCNDSubfolder = 'dataCND/';
dataCNDSubfolder = '';

reRefType = 'Avg'; % or 'Mastoids'
bandpassFilterRange = [0.1,35]; % for MEG
% [1,8] for EEG; 
% Hz (indicate 0 to avoid running the low-pass or high-pass filters or both)
% e.g., [0,8] will apply only a low-pass filter at 8 Hz

downFs = 50;

eegFilenames = dir([dataMainFolder,dataCNDSubfolder,'dataSub*.mat']);
nSubs = length(eegFilenames);

stimIdx = 4; % 1: env; 2: word onset; 3: f0; 4: sp; 5: ap; 6: vuv;
prefix = 'pre_50';

%% Preprocess EEG data
for sub = 1:nSubs
    % Loading EEG data
    eegFilename = [dataMainFolder,dataCNDSubfolder,eegFilenames(sub).name];
    disp(['Loading EEG data: ',eegFilenames(sub).name])
    load(eegFilename,'eeg')
    eeg = cndNewOp(eeg,'Load'); % Saving the processing pipeline in the eeg struct

    % Filtering - LPF (low-pass filter)
    if bandpassFilterRange(2) > 0
        hd = getLPFilt(eeg.fs,bandpassFilterRange(2));
        
        % Filtering each trial/run with a cellfun statement
        eeg.data = cellfun(@(x) filtfilthd(hd,x),eeg.data,'UniformOutput',false);
        
        % Filtering external channels
        if isfield(eeg,'extChan')
            for extIdx = 1:length(eeg.extChan)
                eeg.extChan{extIdx}.data = cellfun(@(x) filtfilthd(hd,x),eeg.extChan{extIdx}.data,'UniformOutput',false);
            end
        end
        
        eeg = cndNewOp(eeg,'LPF');
    end
    
    % Downsampling EEG and external channels
     eeg = cndDownsample(eeg,downFs);
     % stim{1} = cndDownsample(stim{1},downFs);
    
    % Filtering - HPF (high-pass filter)
    if bandpassFilterRange(1) > 0 
        hd = getHPFilt(eeg.fs,bandpassFilterRange(1));
        
        % Filtering EEG data
        eeg.data = cellfun(@(x) filtfilthd(hd,x),eeg.data,'UniformOutput',false);
        
        % Filtering external channels
        if isfield(eeg,'extChan')
            for extIdx = 1:length(eeg.extChan)
                eeg.extChan{extIdx}.data = cellfun(@(x) filtfilthd(hd,x),eeg.extChan{extIdx}.data,'UniformOutput',false);
            end  
        end
        
        eeg = cndNewOp(eeg,'HPF');
    end
    
    % Replacing bad channels
    if isfield(eeg,'chanlocs')
        for tr = 1:length(eeg.data)
            eeg.data{tr} = removeBadChannels(eeg.data{tr}, eeg.chanlocs);
        end
    end
    
    % Re-referencing EEG data
    %eeg = cndReref(eeg,reRefType);
    
    % Removing initial padding (specific to this dataset)
    if isfield(eeg,'paddingStartSample')
        for tr = 1:length(eeg.data)
            eeg.data{tr} = eeg.data{tr}(eeg.paddingStartSample,:);
            for extIdx = 1:length(eeg.extChan)
                eeg.extChan{extIdx}.data = eeg.extChan{extIdx}.data{tr}(1+eeg.paddingStartSample,:);
            end
        end
    end
    
    % Saving preprocessed data
    eegPreFilename = [dataMainFolder,dataCNDSubfolder,prefix,eegFilenames(sub).name];
    disp(['Saving preprocessed EEG data: pre_',eegFilenames(sub).name])
    save(eegPreFilename,'eeg')
end

%% Get average of individual predictions with normal brain data for numSubs subjects.
% Rarely used.

% Loading Stim data

stimIdx = 5;
fileEnd = "_model_32_feature" + stimIdx + ".mat";
modelFilenames = ["sub1" + fileEnd, "sub2" + fileEnd, "sub3" + fileEnd, "sub4" + fileEnd, "sub5" + fileEnd];
avgStruct = struct([]);
for i = 1:20
    wVals = [];
    bVals = [];
    tVals = [];
    
    % Change 4 to 3 and 3 to 2 when the feature is scalar, i.e. stimIdx=3.
    for curr = 1:size(modelFilenames,2)
        currFile = modelFilenames(curr);
        data = load(currFile);
        currModel = data.modelAll(:,i);
        wVals = cat(4,wVals,currModel.w);
        bVals = cat(3, bVals, currModel.b);
        tVals = cat(4,tVals,currModel.t);
    end
    averageModel = struct('w', mean(wVals,4), 'b', mean(bVals,3), 't', mean(tVals,4), 'fs', 64, 'Dir', -1, 'type', 'multi');
    avgStruct = [avgStruct, averageModel];
end

name = "avg" + fileEnd;
modelAll = avgStruct;
save(name, 'modelAll');

%% Get individual prediction with normal brain data for one subject.
% Use for average/individual/rand model.

dirTRF = -1; % Backward TRF model
stimIdx = 3;

% Loading Stim data
stimFilename = 'stim_input_files/meg_audio_files100_4.mat';
load(stimFilename,'stim')

eegPreFilename = "../average_immegSubData.mat";
[stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim, 155);
% For running rand model only
% stimFeature = randomiseObservations(stimFeature, 4); 
% minLen = 2775;
% for x = 1:4
%     stimFeature.data{1,x} = stimFeature.data{1,x}(1:minLen,:);
%     eeg.data{1,x} = eeg.data{1,x}(1:minLen,:);
% end
modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, 0, "avgmeg", 155);
savePred(stimFeature, eeg, modelAll, stimIdx, 0, "avgmeg", 155);

%% MCCA
% For MCCA models.

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 5;
subIdx = 0;

% Loading Stim data
stimFilename = 'meg_audio_files100_4.mat';
% stimFilename = 'stim_input_files/dataStim_32.mat';
load(stimFilename,'stim');

eegPreFilename = 'mcca/subData_meg_imavg.mat';
disp('Loading preprocessed EEG data')

comps = [8,16,32,64,128];
for i = 1:size(comps,2)
    comp = comps(i);
    [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim, comp);
%     stimFeature = randomiseObservations(stimFeature, 4);
    stimFeature.data = flip(stimFeature.data);
    modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, subIdx, "imavg_meg", comp);
    savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, "imavg_meg", comp);
end

%% MCCA individual components

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 3;
subIdx = 0;

for comp = 1:16
    % Loading Stim data
    stimFilename = 'stim_input_files/dataStim_32.mat';
    load(stimFilename,'stim')
    
    eegPreFilename = ['mcca/components/subData_comp',num2str(comp),'.mat'];
    disp('Loading preprocessed EEG data')
    prefix = ['mcca/components/res/comp', num2str(comp)];
    
    [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
    modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, subIdx, prefix);
    savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, prefix);
end

%% MCCA FOR EACH BAND

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 4;
subIdx = 0;

% Loading Stim data
stimFilename = 'stim_input_files/dataStim_32.mat';
load(stimFilename,'stim')

eegPreFilename = 'mcca/subData_pre_128.mat';
disp('Loading preprocessed EEG data')

[stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
thePreds = cell(32,1);
theRs = cell(32,1);
for b = 1:32
    stimFeature2 = stimFeature;
    for i = 1:20
        stimFeature2.data{1,i} = stimFeature2.data{1,i}(:,b);
    end
    modelAll = saveModel(stimFeature2, eeg, dirTRF, stimIdx, subIdx, "");
    [hi1, hi2] = savePred(stimFeature2, eeg, modelAll, stimIdx, subIdx, "");
    theRs{b}= hi1;
    thePreds{b} = hi2;
end


%%

function stimFeature = randomiseObservations(stimFeature, numTrials)
    oldStimFeat = stimFeature;
    r = randperm(numTrials);
    for i = 1:numTrials
        stimFeature.data{i} = oldStimFeat.data{r(i)};
    end
end

function eeg = pickComponents(eeg, components)
    for i = 1:size(eeg.data,2)
        origData = eeg.data{i};
        newData = zeros(size(origData,1), size(components,2));
        for j = 1:size(components,2)
            hello = origData(:,components(j));
            for x = 1:size(hello,1)
                newData(x,j) = hello(x,1);
            end
        end
        eeg.data{i} = newData;
    end

end

function [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim, comps)
    load(eegPreFilename,'eeg'); 

    % only uncomment next line if you want to use a subset of MCCA comps.
    V = uint32(1):uint32(comps);
    eeg = pickComponents(eeg, V);

    % Making sure that stim and neural data have the same length
    stimFeature = stim;
    stimFeature.data = stimFeature.data(stimIdx,:); % chosen stimulus feature
    if eeg.fs ~= stimFeature.fs
        disp('Error: EEG and STIM have different sampling frequency')
        return
    end
    if length(eeg.data) ~= length(stimFeature.data)
        disp('Error: EEG.data and STIM.data have different number of trials')
        return
    end
    for tr = 1:length(stimFeature.data)
        envLen = size(stimFeature.data{tr},1);
        eegLen = size(eeg.data{tr},1);
        minLen = min(envLen,eegLen);
        stimFeature.data{tr} = double(stimFeature.data{tr}(1:minLen,:));
        eeg.data{tr} = double(eeg.data{tr}(1:minLen,:));
    end
    
    % Normalising EEG data
    clear tmpEnv tmpEeg
%     tmpEnv = stimFeature.data{1}; %
    tmpEeg = eeg.data{1};
    for tr = 2:length(stimFeature.data) % getting all values
%         tmpEnv = cat(1,tmpEnv,stimFeature.data{tr}); %
        tmpEeg = cat(1,tmpEeg,eeg.data{tr});
    end
%     normFactorEnv = std(tmpEnv(:)); clear tmpEnv; %
    normFactorEeg = std(tmpEeg(:)); clear tmpEeg;
    for tr = 1:length(stimFeature.data) % normalisation
%         stimFeature.data{tr} = stimFeature.data{tr}/normFactorEnv; %
        eeg.data{tr} = eeg.data{tr}/normFactorEeg;
    end
end


function modelAll = saveModel(stimFeature,eeg, dirTRF, stimIdx, subIdx, prefix, comps)
    clear rAll
    clear modelAll
    
    % TRF - Compute model weights
    disp('Running mTRFcrossval')
    tmin = -200;
    tmax = 600;
    lambdas = [1e-6,1e-4,1e-2,1e0,1e2,1e4];
%     tmin = -100;
%     tmax = 400;
%     lambdas = [1e0,1e2,1e4]; 

    for obsIdx = 1:length(stimFeature.data)
        trainStim = stimFeature.data;
        trainStim(:,obsIdx) = [];
        trainResp = eeg.data;
        trainResp(:,obsIdx) = [];
        
        [stats,t] = mTRFcrossval(trainStim,trainResp,eeg.fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
        [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));
        bestLambda
        disp(['r = ',num2str(maxR)])
        %bestLambda=6;
    
        disp('Running mTRFtrain')
        model = mTRFtrain(trainStim,trainResp,eeg.fs,dirTRF,tmin,tmax,lambdas(bestLambda),'verbose',0);
        modelAll(obsIdx) = model;
    end
    name = prefix + "sub" + subIdx + "_model_100_" + num2str(comps) + "_feature" + stimIdx + ".mat";
    save(name,"modelAll");
end


function [rVals, predAll] = savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, prefix, comps)
    clear predAll
    rVals = [];
    numObs = length(stimFeature.data);
    for obsIdx = 1:numObs
        testStim = stimFeature.data(:,obsIdx);
        testResp = eeg.data(:,obsIdx);
        model = modelAll(obsIdx);
        [pred,test] = mTRFpredict(testStim,testResp,model,'verbose',0);
    
        % Plot test accuracy
        rVals = cat(1,rVals,test.r);
        subplot(2,2,4), bar(1,0.07006), hold on, bar(2,test.r), hold off
        set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
        title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
        predAll{obsIdx} = pred;
    end

    name = prefix + "sub" + subIdx + "_pred_100_" + num2str(comps) + "_feature" + stimIdx + ".mat";
    save(name,"predAll");
    name = prefix + "sub" + subIdx + "_rvals_100_" + num2str(comps) + "_feature" + stimIdx + ".mat";
    save(name,"rVals");
end


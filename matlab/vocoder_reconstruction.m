% Adapted from CNSP-Workshop 2021
% https://cnsp-workshop.github.io/website/index.html
% Author: Giovanni M. Di Liberto
% Copyright 2021 - Giovanni Di Liberto
%                  Nathaniel Zuk
%                  Michael Crosse
%                  Aaron Nidiffer
%                  (see license file for details)
% Last update: 11 February 2022

clear all
close all

addpath libs/cnsp_utils
addpath libs/cnsp_utils/cnd
addpath libs/mTRF-Toolbox_v2/mtrf
addpath libs/NoiseTools
addpath eeglab2021.1


%% Parameters - Natural speech listening experiment
dataMainFolder = '../datasets/LalorNatSpeech/';
% dataMainFolder = '../datasets/LalorNatSpeechReverse/';
dataCNDSubfolder = 'dataCND/';

reRefType = 'Avg'; % or 'Mastoids'
bandpassFilterRange = [1,8]; % Hz (indicate 0 to avoid running the low-pass
                          % or high-pass filters or both)
                          % e.g., [0,8] will apply only a low-pass filter
                          % at 8 Hz
% changed by sophie
downFs = 32; % Hz. *** fs/downFs must be an integer value ***

eegFilenames = dir([dataMainFolder,dataCNDSubfolder,'dataSub*.mat']);
nSubs = length(eegFilenames);

stimIdx = 3; % 1: env; 2: word onset; 3: f0; 4: sp; 5: ap; 6: vuv;

%% Preprocess EEG - Natural speech listening experiment
for sub = 1:nSubs
    % Loading EEG data
    eegFilename = [dataMainFolder,dataCNDSubfolder,eegFilenames(sub).name];
    disp(['Loading EEG data: ',eegFilenames(sub).name])
    load(eegFilename,'eeg')
    eeg = cndNewOp(eeg,'Load'); % Saving the processing pipeline in the eeg struct

    % Filtering - LPF (low-pass filter)
    if bandpassFilterRange(2) > 0
        hd = getLPFilt(eeg.fs,bandpassFilterRange(2));
        
        % Filtering each trial/run with a for loop
%         for ii = 1:length(eeg.data)
%             eeg.data{ii} = filtfilthd(hd,eeg.data{ii});
%         end
        
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
    eeg = cndReref(eeg,reRefType);
    
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
    eegPreFilename = [dataMainFolder,dataCNDSubfolder,'pre_',eegFilenames(sub).name];
    disp(['Saving preprocessed EEG data: pre_',eegFilenames(sub).name])
    save(eegPreFilename,'eeg')
end

%% 
% Get individual predictions with normal eeg for numSubs subjects.

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 5;
numSubs = 5;

% Loading Stim data
stimFilename = 'dataStim_32.mat';
disp(['Loading stimulus data: ','dataStim_32.mat'])
load(stimFilename,'stim')

% Get models for: 1 feature, [subIdx] subjects, all observations.
for subIdx = 1:numSubs
    % Loading preprocessed EEG
    eegPreFilename = [dataMainFolder,dataCNDSubfolder,'pre_',eegFilenames(subIdx).name];
    disp(['Loading preprocessed EEG data: pre_',eegFilenames(subIdx).name])
    [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
    modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, subIdx, "");
    savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, "");
end

%%
%% 
% normal individual.

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 5;

% Loading Stim data
stimFilename = 'dataStim_32.mat';
disp(['Loading stimulus data: ','dataStim_32.mat'])
load(stimFilename,'stim')

% Loading preprocessed EEG
eegPreFilename = 'avgSubData';
disp(['Loading preprocessed EEG data: ',eegPreFilename])
[stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, 0, "avg2");
savePred(stimFeature, eeg, modelAll, stimIdx, 0, "avg2");

%%
% Get average of individual predictions with normal eeg for numSubs subjects.

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

%%

% Get individual prediction with normal eeg for one subject.

stimIdx = 3;

% Loading Stim data
stimFilename = 'dataStim_32.mat';
disp(['Loading stimulus data: ','dataStim_32.mat'])
load(stimFilename,'stim')

eegPreFilename = [dataMainFolder,dataCNDSubfolder,'pre_',eegFilenames(1).name];
disp(['Loading preprocessed EEG data: pre_',eegFilenames(1).name])
[stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
savePred(stimFeature, eeg, modelAll, stimIdx, 1, "avg_");

%% 
%  MCCA

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 5;
subIdx = 0;

% Loading Stim data
stimFilename = 'stim_input_files/dataStim_32.mat';
disp(['Loading stimulus data: ','dataStim_32.mat'])
load(stimFilename,'stim')

eegPreFilename = 'mcca/subData_pre_128.mat';
disp('Loading preprocessed EEG data')

[stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, subIdx, "cmb");
savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, "cmb");

%% testing an avg thing

stimIdx = 3;

% Loading avg model
avgModelFilename = ['avg_sub_results/avg_model_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','avg_model_32_feature', num2str(stimIdx), '.mat']);
avgModel = load(avgModelFilename,'modelAll').modelAll;

curr = zeros(20,19);
for sub = 1:19
    clear stim stimFeature eeg
    % Loading Stim data
    stimFilename = 'dataStim_32.mat';
    disp(['Loading stimulus data: ','dataStim_32.mat'])
    load(stimFilename,'stim')
    
    % eeg
    eegPreFilename = ['../datasets/LalorNatSpeech/dataCND/pre_dataSub',num2str(sub),'.mat'];
    disp('Loading preprocessed EEG data')
    [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
    currVals = savePred(stimFeature, eeg, avgModel, stimIdx, sub, "");
    curr(:,sub) = currVals;
end

rVals = zeros(20,1);
for i=1:20
    rVals(i,:) = mean(curr(i,:));
end

name =  "avg_rvals_32_feature" + stimIdx + ".mat";
save(name,"rVals");

%% 
%  MCCA individual components

clear modelAll
dirTRF = -1; % Backward TRF model
stimIdx = 3;
subIdx = 0;

for comp = 1:16
    % Loading Stim data
    stimFilename = 'stim_input_files/dataStim_32.mat';
    disp(['Loading stimulus data: ','dataStim_32.mat'])
    load(stimFilename,'stim')
    
    eegPreFilename = ['mcca/comps/subData_comp',num2str(comp),'.mat'];
    disp('Loading preprocessed EEG data')
    prefix = ['mcca/comps/res/comp', num2str(comp)];
    
    [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim);
    modelAll = saveModel(stimFeature, eeg, dirTRF, stimIdx, subIdx, prefix);
    savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, prefix);
end

%% testing thing

    stimIdx = 3;
    eegPreFilename = 'mcca/subData.mat';    
    load(eegPreFilename,'eeg');
    eeg = pickComponents(eeg, [1,3,5,7,8,9,10,15,16]);
    size(eeg.data{1})
    
%%

function eeg = pickComponents(eeg, components)
    for i = 1:20
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

function [stimFeature, eeg] = preprocessData(eegPreFilename, stimIdx, stim)
    load(eegPreFilename,'eeg');

    % only uncomment next line if you want to use a subset of MCCA comps.
    V = uint32(1):uint32(64);
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
    tmpEnv = stimFeature.data{1};
    tmpEeg = eeg.data{1};
    for tr = 2:length(stimFeature.data) % getting all values
        tmpEnv = cat(1,tmpEnv,stimFeature.data{tr});
        tmpEeg = cat(1,tmpEeg,eeg.data{tr});
    end
    normFactorEnv = std(tmpEnv(:)); clear tmpEnv;
    normFactorEeg = std(tmpEeg(:)); clear tmpEeg;
    for tr = 1:length(stimFeature.data) % normalisation
        stimFeature.data{tr} = stimFeature.data{tr}/normFactorEnv;
        eeg.data{tr} = eeg.data{tr}/normFactorEeg;
    end
end


function modelAll = saveModel(stimFeature,eeg, dirTRF, stimIdx, subIdx, prefix)
    clear rAll
    clear modelAll
    
    % TRF - Compute model weights
    disp('Running mTRFcrossval')
    tmin = -200;
    tmax = 600;
    lambdas = [1e-6,1e-4,1e-2,1e0,1e2,1e4];

    for obsIdx = 1:length(stimFeature.data)
        trainStim = stimFeature.data;
        trainStim(:,obsIdx) = [];
        trainResp = eeg.data;
        trainResp(:,obsIdx) = [];
        
        [stats,t] = mTRFcrossval(trainStim,trainResp,eeg.fs,dirTRF,tmin,tmax,lambdas,'verbose',0);
        [maxR,bestLambda] = max(squeeze(mean(mean(stats.r,1),3)));
        disp(['r = ',num2str(maxR)])
    
        disp('Running mTRFtrain')
        model = mTRFtrain(trainStim,trainResp,eeg.fs,dirTRF,tmin,tmax,lambdas(bestLambda),'verbose',0);
        modelAll(obsIdx) = model;
    end
    name = prefix + "sub" + subIdx + "_model_32_feature" + stimIdx + ".mat";
    save(name,"modelAll");
end


function rVals = savePred(stimFeature, eeg, modelAll, stimIdx, subIdx, prefix)
    clear predAll
    rVals = [];
    for obsIdx = 1:length(stimFeature.data)
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

    name = prefix + "sub" + subIdx + "_pred_32_feature" + stimIdx + ".mat";
    save(name,"predAll");
    name = prefix + "sub" + subIdx + "_rvals_32_feature" + stimIdx + ".mat";
    save(name,"rVals");
end


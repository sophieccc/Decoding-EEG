clear variables
close all

addpath libs/cnsp_utils
addpath libs/cnsp_utils/cnd
addpath libs/mTRF-Toolbox_v2/mtrf
addpath libs/NoiseTools
addpath eeglab2021.1

%%
stimIdx = 3;

%% Compare r values of bands


mcca_r_file = ['results_mcca/cmb/cmb_rvals_32_feature4.mat'];
mccaR = load(mcca_r_file,'rVals').rVals;
curr = zeros(32,1);

for band = 1:32
    theMean = mean(mccaR(:,band));
    curr(band,1) = theMean;
end

h = bar(curr);
title('Mean R Value per Spectrogram Band'), xlabel('Band'), ylabel('R')
run prepExport.m
%%
% Loading mcca prediction data
mccaFile = ['cmbsub0_pred_32_feature3.mat'];
pred = load(mccaFile,'predAll').predAll{1,1};

% Loading stimulus data (orig)
origFile = 'stim_input_files/audio_files.mat';
orig = load(origFile,'stim').stim;
origFeat = orig.data{3,1};

% Loading stimulus data (with filt)
stimFile = 'stim_input_files/fixedStim32.mat';
stim = load(stimFile,'stim').stim;
feature = stim.data{3,1};
feature = feature(1:size(pred,1),:);

% Loading stimulus data (without filt)
stimFile2 = 'newStim32.mat';
stim2 = load(stimFile2,'stim').stim;
feature2 = stim2.data{3,1};
feature2 = feature2(1:size(pred,1),:);

% Loading stimulus data (before all this)
stimFile3 = 'stim_input_files/dataStim_32.mat';
stim3 = load(stimFile3,'stim').stim;
feature3 = stim3.data{3,1};
feature3 = feature3(1:size(pred,1),:);
size(feature3)

subplot(2,1,1);
plot(pred);
subplot(2,1,2); 
plot(feature3);
% subplot(3,1,3); 
% plot(feature3);

%%
% Compare no. of compnents

% feature 3: 
 x = [1,4,8, 16,32,128];
  y = [0.1629,0.18137,0.18514, 0.18532,0.18357,0.17817];
% before fixing f0: y = [0.17051,0.18946,0.1937,0.19223,0.1876];

% feature 4: 
% x = [1,3,7,16,128];
% y = [0.11822,0.12549,0.12734,0.12375,0.10694];

% feature 5: 
%x = [1,7,16,32,64,128];
%y = [0.15635,0.19133,0.19229,0.19209, 0.19627,0.19477];

plot(x,y,'linewidth',2),
title('R Values for MCCA Component Combinations [Feature 3]')
xlabel('No. of Components')
ylabel('R')
axis square, grid on
run prepExport.m

%%
r_file = 'results_meg/megsub5_rvals_100_32_feature5.mat'; 
rVals = load(r_file,'rVals').rVals;

curr = zeros(20,1);
for row = 1:20
    total = 0;
    for col = 1:size(rVals,2)
        val = rVals(row,col);
        total = total + abs(val);
    end
    total = total / size(rVals,2);
    curr(row, 1) = total;
end

h = bar(curr);
run prepExport.m

%testAccuracy(curr(:, model));
disp(['Mean r = ',num2str(mean(curr(:,1)))])

%% Compare r values

stimIdx = 3;

indiv_r_file = ['results_indiv/sub1_rvals_32_feature', num2str(stimIdx), '.mat'];
indivR = load(indiv_r_file,'rVals').rVals;

avg_r_file = ['results_avg/avg2sub0_rvals_32_feature', num2str(stimIdx), '.mat'];
avgR = load(avg_r_file,'rVals').rVals;

mcca_r_file = ['results_mcca/cmb/cmb_rvals_32_feature', num2str(stimIdx), '.mat'];
mccaR = load(mcca_r_file,'rVals').rVals;

data = {indivR;avgR;mccaR};
curr = zeros(20,3);

for i = 1:3
    currData = data{i,:};
    for row = 1:20
        total = 0;
        for col = 1:size(currData,2)
            val = currData(row,col);
            total = total + abs(val);
        end
        total = total / size(currData,2);
        curr(row, i) = total;
    end
end

h = bar(curr);
set(h, {'DisplayName'}, {'Individual','Average','MCCA'}')
title('Mean R Value per Observation for F0'), xlabel('Observation'), ylabel('R')
run prepExport.m
legend()

names = ["Individual","Average", "MCCA"];
for model = 1:3
    %testAccuracy(curr(:, model));
    msg = "Mean r for " + names(model) + " [feature " + num2str(stimIdx) +"] = " + num2str(mean(curr(:,model)));
    disp(msg)
end

%% Compare r values overall
means = [];
for idx = 3:5
    rnd_r_file = ['results_mcca/random/rndsub0_rvals_32_feature', num2str(idx), '.mat'];
    rndR = load(rnd_r_file,'rVals').rVals;

    indiv_r_file = ['results_indiv/sub1_rvals_32_feature', num2str(idx), '.mat'];
    indivR = load(indiv_r_file,'rVals').rVals;
    
    avg_r_file = ['results_avg/avg2sub0_rvals_32_feature', num2str(idx), '.mat'];
    avgR = load(avg_r_file,'rVals').rVals;
    
    mcca_r_file = ['results_mcca/cmb/cmb_rvals_32_feature', num2str(idx), '.mat'];
    mccaR = load(mcca_r_file,'rVals').rVals;
    

    data = {rndR;indivR;avgR;mccaR};
    curr = zeros(20,4);

    for i = 1:4
        currData = data{i,:};
        for row = 1:20
            total = 0;
            for col = 1:size(currData,2)
                val = currData(row,col);
                total = total + abs(val);
            end
            total = total / size(currData,2);
            curr(row, i) = total;
        end
    end
    
    means = [means; [mean(curr(:,1)) mean(curr(:,2)) mean(curr(:,3)) mean(curr(:,4))]];
end

x = categorical(["F0" "Spectrogram" "Aperiodicity"]);
h = bar(x,means);
set(h, {'DisplayName'}, {'Random Model','Individual Model','Average Model','MCCA Model'}')
title('Mean R Value per Stimulus Feature'), xlabel('Stimulus Feature'), ylabel('R')
run prepExport.m
legend() 

%% evaluating individual components

names = ["F0", "Spectrogram", "Aperiodicity"];
for stimIdx = 3:5
    subplot(3,1,stimIdx-2);
    curr = zeros(20,16);
    
    for i = 1:16
        filename = ['mcca/components/res/comp', num2str(i), 'sub0_rvals_32_feature', num2str(stimIdx),'.mat'];
        rVals = load(filename,'rVals').rVals;
        currData = rVals;
        for row = 1:20
            total = 0;
            for col = 1:size(currData,2)
                val = currData(row,col);
                total = total + abs(val);
            end
            total = total / size(currData,2);
            curr(row, i) = total;
        end
    end
    
    means = zeros(16,1);
    for j = 1:16
        %testAccuracy(curr(:, model));
        disp(['Mean r for component ', num2str(j), ' [feature ', num2str(stimIdx), '] = ',num2str(mean(curr(:,j)))])
        means(j,1) = mean(curr(:,j));
    end
    
    h = bar(means);
    theTitle = "";
    if stimIdx == 3
        theTitle = 'Mean R Values Per Component';
    end
    
    title(theTitle), xlabel(names(stimIdx-2)+ ' Component'), ylabel('R')
    run prepExport.m
end

%%

% Loading mcca model data
modelFilename = ['results_mcca/cmb/cmb_model_32_feature', num2str(stimIdx), '.mat'];
mccaModel = load(modelFilename,'modelAll').modelAll;

% Loading avg model data
avgModelFilename = ['results_avg/avg2sub0_model_32_feature', num2str(stimIdx), '.mat'];
avgModel = load(avgModelFilename,'modelAll').modelAll;


mccaAvgModel = plotTRF(mccaModel);
plotGFP(mccaAvgModel);

%avgAvgModel = plotTRF(avgModel);
%plotGFP(avgAvgModel);

%%

% Loading mcca prediction data
mccaFile = ['results_mcca/cmb/cmb_pred_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','sub0_pred_32_feature', num2str(stimIdx), '.mat']);

% Loading avg prediction data
avgFile = ['results_avg/avg2sub0_pred_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','avg2sub0_pred_32_feature', num2str(stimIdx), '.mat']);

% Loading stimulus data
stimFile = 'stim_input_files/dataStim_32.mat';
disp('Loading model data: dataStim_32.mat');
stim = load(stimFile,'stim').stim;

valsMCCA = [];
valsAvg = [];
for i = 1:20
    mccaPred = load(mccaFile,'predAll').predAll{i};

    avgPred = load(avgFile,'predAll').predAll{i};
    avgPred = avgPred(1:size(mccaPred,1),:);

    feature = stim.data{stimIdx,i};
    feature = feature(1:size(mccaPred,1),:);
    
    mccaCorr = corrcoef(feature,mccaPred);
    mccaCorr
    otherCorr = mTRFevaluate(feature,mccaPred);
    mean(otherCorr)
    valsMCCA = [valsMCCA mccaCorr(1,2)];
    avgCorr = corrcoef(feature,avgPred);
    valsAvg = [valsAvg avgCorr(1,2)];
end

disp(['Mean corr for MCCA = ',num2str(mean(valsMCCA))])
disp(['Mean corr for Avg = ',num2str(mean(valsAvg))])

%% Comparing f0 correlation 

% Loading mcca prediction data
mccaFile = ['nofiltnorm_pred_32_feature3.mat'];

% Loading stimulus data
stimFile = 'newStim32.mat';
stim = load(stimFile,'stim').stim;

valsMCCA = [];
oldValsMCCA = [];
for i = 1:20
    mccaPred = load(mccaFile,'predAll').predAll{i};
    feature = stim.data{stimIdx,i};
    feature = feature(1:size(mccaPred,1),:);
    
    newPred = [];
    newFeat = [];
    theNum = quantile(feature,0.5);
    for j = 1:size(mccaPred,1)
        if feature(j) >= theNum
            newPred = [newPred mccaPred(j)];
            newFeat = [newFeat feature(j)];
        end
    end

    oldMccaCorr = corrcoef(feature,mccaPred);
    oldValsMCCA = [oldValsMCCA oldMccaCorr(1,2)];

    mccaCorr = corrcoef(newFeat,newPred);
    valsMCCA = [valsMCCA mccaCorr(1,2)];
end

disp(['Mean corr for MCCA = ',num2str(mean(valsMCCA))])
disp(['Old mean corr for MCCA = ',num2str(mean(oldValsMCCA))])

%% FIND BEST AUDIO

r_file = 'results_mcca/nonorm/nonorm_rvals_32_feature4.mat'; 
rVals = load(r_file,'rVals').rVals;

curr = zeros(20,1);
for row = 1:20
    total = 0;
    for col = 1:size(rVals,2)
        val = rVals(row,col);
        total = total + abs(val);
    end
    total = total / size(rVals,2);
    curr(row, 1) = total;
end

curr

%% FIND BEST PART OF AUDIO

stimIdx = 3;
% Loading mcca prediction data
mccaFile = ['results_mcca/nonorm/nonorm_pred_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','nonorm_pred_32_feature', num2str(stimIdx), '.mat']);

% Loading stimulus data
stimFile = 'stim_input_files/dataStim_32.mat';
disp('Loading model data: dataStim_32.mat');
stim = load(stimFile,'stim').stim;

valsMCCA = [];
a = [1:12];
b = a.*400;
b = [1 b];
for j = 1:size(b,2)
    firstNum = b(j);
    secondNum = firstNum + 399;
    mccaPred = load(mccaFile,'predAll').predAll{9};

    feature = stim.data{stimIdx,9};
    feature = feature(1:size(mccaPred,1),:);
    mccaPred = mccaPred(firstNum:secondNum);
    feature = feature(firstNum:secondNum);
    
    mccaCorr = corrcoef(feature,mccaPred);
    valsMCCA = [valsMCCA mccaCorr(1,2)];
end

%%
stimIdx = 3;
% Loading mcca prediction data
mccaFile = ['results_mcca/nonorm/nonorm_pred_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','nonorm_pred_32_feature', num2str(stimIdx), '.mat']);

% Loading stimulus data
stimFile = 'stim_input_files/dataStim_32.mat';
disp('Loading model data: dataStim_32.mat');
stim = load(stimFile,'stim').stim;

 mccaPred = load(mccaFile,'predAll').predAll{9};
 mccaPred = mccaPred(4400:4499);
 feature = feature(1:size(mccaPred,1),:);
 feature = feature(firstNum:secondNum);

%%

function testAccuracy(rAll)
    for i = 1:20
        r = rAll(i);
        subplot(2,2,4), bar(1,0.07006), hold on, bar(2,r), hold off
        set(gca,'xtick',1:2,'xticklabel',{'Val.','Test'}), axis square, grid on
        title('Model Performance'), xlabel('Dataset'), ylabel('Correlation')
    end
    drawnow;
end

function plotGFP(avgModel)
    % Plot GFP
    subplot(1,2,2)
    tmin = -200;
    tmax = 600;
    mTRFplot(avgModel,'gfp',[],'all',[tmin+50,tmax-50]);
    title('Global Field Power')
    run prepExport.m
    
    drawnow;
end

function [avgModel] = plotTRF(modelAll)
    % todo: don't hardcode (what do you mean? where would I get it?)
    tmin = -200;
    tmax = 600;

    % Plot average TRF
    avgModel = mTRFmodelAvg(modelAll,0);

    % Plot avg TRF model
    subplot(1,2,1)

    el = 1;
    mTRFplot(avgModel,'trf',[],el, [tmin+50,tmax-50]);

    %plot(avgModel.t,squeeze(avgModel.w))
    title('Envelope avgTRF')
    ylabel('Magnitude (a.u.)')

    axis square
    run prepExport.m
    drawnow;
end

function thisSTI = calcSTI(filename)
% Gets STI given audio filename    
    info = audioinfo(filename);
    fs = info.SampleRate;
    data = audioread(filename);
    thisSTI = STI(data, fs);
    hold on
    plot(data, 'DisplayName', filename)
end

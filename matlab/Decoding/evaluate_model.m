% Many methods of evaluating our results.

clear variables
close all

addpath libs/cnsp_utils
addpath libs/cnsp_utils/cnd
addpath libs/mTRF-Toolbox_v2/mtrf
addpath libs/NoiseTools
addpath eeglab2021.1

%%
stimIdx = 3;

%%
imagined=load('stim_output_files/meg/immcca_ld.mat').stim.data(4,:);
listened=load('stim_output_files/meg/mcca_ld.mat').stim.data(4,:);
standard=load('stim_output_files/meg/theStandard.mat').stim.data(4,:);


[h,p] = ttest(imagined{1,1},listened{1,2});
mean(p)

[h,p] = ttest(imagined{1,3},listened{1,4});
mean(p)

[h,p] = ttest(imagined{1,1},listened{1,3});
mean(p)

[h,p] = ttest(imagined{1,1},listened{1,4});
mean(p)

[h,p] = ttest(imagined{1,2},listened{1,3});
mean(p)

[h,p] = ttest(imagined{1,2},listened{1,4});
mean(p)

%% Values gotten from variations of code section above.

mean([0.0806 0.2452]) % p values for same listened audios [1and2 3and4]
mean([0.0178 0.0013 0.00000051899 0.00011318]) % p values for diff listened audios
mean([0.3186 0.0432]) % p values for same imagined audios
mean([0.0281 0.0029 0.00000054234 0.00061125]) % p values for diff imagined audios

mean([0.1340 0.0295 0.0380 0.2225]) % p values for same listened + imagined 
mean([0.0028 0.0237 0.0263 0.0301 0.0184 0.0427 0.000073082 0.0040]) % p values for diff listened + imagined 

%% Creating audio metrics bar plot

% listened
% data = [
%     normalize([0.2860018121365055 0.29875158413694586 0.30163349930444106 0.2962736720206629 0.302096653804104 0.3363945585658703 0.3244315372119114], 'norm', 1); 
%     normalize([-0.49478337024439634 -0.42699465833957445 -0.3952664691058791 -0.24926329903408456 -0.7434235445197486 -0.2505483029656429 -0.572950566410872], 'norm', 1); 
%     normalize([9.253541668180528 9.271305123569912 9.24131583236917 9.185993289866639 9.245784063246063 9.241772873971428 9.260791671982108], 'norm',1)
%     ];

% imagined
% data = [
%     normalize([0.28311348903938943 0.27425094406127815 0.2604431636850755 0.2873664848765392 0.2953004393601648 0.2839023532185543 0.2704607289362801], 'norm', 1); 
%     normalize([-0.43118917874349283 -0.5157774509241625 -0.502065976333567 -0.2843252027044656 -0.8242783437791662 -0.29205865855987667 -0.6051727327093097], 'norm', 1); 
%     normalize([9.192779959450256 9.222072708501592 9.286020376565073 9.17191674467727 9.253252379373702 9.368284566368505 9.388788713545551], 'norm', 1)
%     ];

% EEG
% data = [
%     normalize([0.38096 0.37229 0.37353 0.39105 0.40497], 'norm', 1); 
%     normalize([0.39607 0.1967 0.22224 0.18571 0.6048], 'norm', 1); 
%     normalize([7.72078 7.79697 7.84664 7.7501 7.78162], 'norm', 1)
%     ];

% improvements on imagined MEG
data = [
    normalize([0.28737 0.33490036982735816 0.6509857052166141 0.2419145901307067 0.2953776591330873 0.312810 0.47844579581756197 0.4774], 'norm', 1); 
    normalize([-0.28433 -0.5994865745717542 4.510735854746559 -0.20559755812410577 -0.8256425575875784 0.59812 0.390084213810745 0.88643], 'norm', 1); 
    normalize([9.17192 9.244102170793273 3.7122182672826556 9.21233742972544 9.253478084060244 7.88534 7.9016523406063275 7.42453], 'norm', 1)
    ];

x = categorical([ "STOI" "fwSNRSeg" "CD"]);

h = bar(x,data);
hold on

title('Audio Metric Results Across Imagined MEG Models'), xlabel('Audio Metric'), ylabel('Value')
run prepExport.m
% MEG
% legend('Random Model','Individual Model','Average Model','MCCA Model', 'MCCA Model (Norm F0)', 'MCCA Model (Small λ)', 'MCCA Model (Small λ, Norm F0)') 
% EEG
% legend('Random Model','Individual Model','Average Model', 'Normalised MCCA Model', 'MCCA Model');

%Improvements
legend('MCCA','Real F0','Real Spectrogram','Real Aperiodicity', 'Average Trials', 'Scale Features', 'Modulate Envelope', 'Scale + Envelope') 
%% Looking at correlations between MCCA components and EEG electrodes.

mcca = load('mcca/subData_pre_128.mat','eeg').eeg; 
mccaData = mcca.data;
theSize = size(mccaData{1,1},1);
corrMat  = zeros(128,128);
for sub = 1:19
    fileName = ['../datasets/LalorNatSpeech/dataCND/pre_dataSub',num2str(sub),'.mat'];
    eeg = load(fileName,'eeg').eeg;
    eegData = eeg.data;
    for trial = 1:20
        for comp = 1:128
            for elec = 1:128
                    mccaNums = mccaData{1,trial}(comp,:);
                    eegNums = eegData{1,trial}(elec,:);
                    r = corrcoef(mccaNums,eegNums);
                    corrMat(comp,elec) = corrMat(comp,elec) + abs(r(1,2));
            end
        end
    end
end
corrMat = corrMat/(19 * 20);

compToElec = zeros(128,2);
elecToComp = zeros(128,2);
for i = 1:128
    [compToElec(i,1), compToElec(i,2)] = max(corrMat(i,:));
    [elecToComp(i,1), elecToComp(i,2)] = max(corrMat(:,i));
end

[n, bin] = hist(compToElec(:,2), unique(compToElec(:,2)));
[hype,idx] = sort(-n);
hi = bin(idx);

imagesc(corrMat);

[n2, bin2] = hist(elecToComp(:,2), unique(elecToComp(:,2)));
[hype2,idx2] = sort(-n2);
hi2 = bin2(idx2);
%% paired t-test avg vs mcca

num_trials = 4;
comps = [8,16,8];
for idx = 3:5
    comp = num2str(comps(idx-2));
    avg_r_file = ['results_meg/avg_mod/average_immegsub5_rvals_100_155_feature', num2str(idx), '.mat'];
    avgR = load(avg_r_file,'rVals').rVals;
    
    mcca_r_file = ['results_meg/not_ld/imavg_megsub0_rvals_100_',comp,'_feature', num2str(idx), '.mat'];
%     if idx == 5
%         mcca_r_file = ['results_mcca/cmb/cmb_rvals_32_feature', num2str(idx), '.mat'];
%     end
    mccaR = load(mcca_r_file,'rVals').rVals;

    data = {avgR;mccaR};
    curr = zeros(num_trials,2);

    for i = 1:2
        currData = data{i,:};
        for row = 1:num_trials
            total = 0;
            for col = 1:size(currData,2)
                val = currData(row,col);
                total = total + abs(val);
            end
            total = total / size(currData,2);
            curr(row, i) = total;
        end
    end
    
    [h,p] = ttest(curr(:,1),curr(:,2));
    disp(p)
    disp(h)
end

%% stimulus vs reconstruction 

stimFile = 'stim_input_files/dataStim_32.mat';
disp('Loading model data: dataStim_32.mat');
stim = load(stimFile,'stim').stim;
figure;

% feature 3:
tiledlayout(2,1)
nexttile
stimData = stim.data{3,1};
plot(stimData(1000:1299)),


title('Stimulus F0')
xlabel('Time (Samples)')
ylabel('F0 (Hz)')
grid on
run prepExport.m

nexttile
mcca_file = ['results_mcca/nonorm/nonorm_pred_32_feature3.mat'];
mccaPred = load(mcca_file,'predAll').predAll{1,1};
plot(mccaPred(1000:1299)),

title('Reconstructed F0')
xlabel('Time (Samples)')
ylabel('F0 (Hz)')
grid on
run prepExport.m


figure;
% feature 4: 
tiledlayout(2,1)
nexttile
stimData = stim.data{4,1};
imagesc(stimData(1000:1099,:)'),

title('Stimulus Spectrogram')
xlabel('Time (Samples)')
ylabel('Frequency Band')

nexttile
mcca_file = ['results_mcca/nonorm/nonorm_pred_32_feature4.mat'];
mccaPred = load(mcca_file,'predAll').predAll{1,1};
mccaPred(mccaPred<0) = 0;
imagesc(mccaPred(1000:1099,:)'),

title('Reconstructed Spectrogram')
xlabel('Time (Samples)')
ylabel('Frequency Band')


figure;
% feature 5: 
tiledlayout(2,1)
nexttile
stimData = stim.data{5,1};
imagesc(stimData(1000:1099,:)'),

title('Stimulus Aperiodicity')
xlabel('Time (Samples)')
ylabel('Frequency Band')

nexttile
mcca_file = ['results_mcca/nonorm/nonorm_pred_32_feature5.mat'];
mccaPred = load(mcca_file,'predAll').predAll{1,1};
mccaPred(mccaPred<0) = 0;
imagesc(mccaPred(1000:1099,:)'),

title('Reconstructed Aperiodicity')
xlabel('Time (Samples)')
ylabel('Frequency Band')

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
mccaFile = ['results_mcca/cmb/cmb_pred_32_feature3.mat'];
pred = load(mccaFile,'predAll').predAll{1,1};

% Loading stimulus data (orig)
origFile = 'stim_input_files/audio_files.mat';
orig = load(origFile,'stim').stim;
origFeat = orig.data{3,1};

% Loading stimulus data (with filt)
stimFile = 'stim_input_files/filteredStim_32.mat';
stim = load(stimFile,'stim').stim;
feature = stim.data{3,1};
feature = feature(1:size(pred,1),:);

% Loading stimulus data (without filt)
stimFile2 = 'stim_input_files/dataStim_32.mat';
stim2 = load(stimFile2,'stim').stim;
feature2 = stim2.data{3,1};
feature2 = feature2(1:size(pred,1),:);

subplot(2,1,1);
plot(pred);
subplot(2,1,2); 
plot(feature3);
% subplot(3,1,3); 
% plot(feature3);

%%
% Compare no. of components

% feature 3: 
% x = [1,4,8, 16,32,128];
% %y = [0.1629,0.18137,0.18514, 0.18532,0.18357,0.17817];
% y = [0.17051,0.18946, 0.1925, 0.1937,0.19223,0.1876]; % before fixing f0:
subplot(1,3,1);
plot(x,y,'linewidth',2),
title('R Values for F0')
xlabel('No. of Components')
ylabel('R')
axis square, grid on
run prepExport.m

% feature 4: 
% x1 = [1,3,7,16,128];
% y1 = [0.11822,0.12549,0.12734,0.12375,0.10694];
subplot(1,3,2);
plot(x1,y1,'linewidth',2),
title('R Values for Spectrogram')
xlabel('No. of Components')
ylabel('R')
axis square, grid on
run prepExport.m

% feature 5: 
% x2 = [1,7,16,32,64,128];
% y2 = [0.15635,0.19133,0.19229,0.19209, 0.19627,0.19477];
subplot(1,3,3);
plot(x2,y2,'linewidth',2),
title('R Values for Aperiodicity')
xlabel('No. of Components')
ylabel('R')
axis square, grid on
run prepExport.m

%%
% EEG
% x = [1,4,8, 16,32,128];
% y = [0.17051,0.18946, 0.1925, 0.1937,0.19223,0.1876];
% x1 = [1,3,7,16,128];
% y1 = [0.11822,0.12549,0.12734,0.12375,0.10694];
% x2 = [1,7,16,32,64,128];
% y2 = [0.15635,0.19133,0.19229,0.19209, 0.19627,0.19477];

% MEG Listening
x = [8,16,32,64,128];
y = [0.45404,0.47991, 0.49717, 0.51098,.51625];
x1 = [8,16,32,64,128];
y1 = [0.16518,0.18473,0.1773,0.18469,0.18771];
x2 = [8,16,32,64,128];
y2 = [0.44645,0.46883,0.47879,0.48833, 0.48918];

% MEG Imagining
% x = [8,16,32,64,128];
% y = [0.16543,0.20917, 0.34873, 0.33389,0.32932];
% x1 = [8,16,32,64,128];
% y1 = [0.095811,0.091553,0.11052,0.097777,0.10461];
% x2 = [8,16,32,64,128];
% y2 = [0.17479,0.20132,0.33959,0.32567, 0.31982];

plot(x,y,'linewidth',2),
hold on 
plot(x1,y1,'linewidth',2),
hold on 
plot(x2,y2,'linewidth',2),
hold off 
title('R Values for Combinations of Components [Listening]')
xlabel('No. of Components')
ylabel('R')
axis square, grid on
run prepExport.m
legend('F0','Spectrogram','Aperiodicity')

%%

comps = [128];
for i = 1:size(comps,2)
    comp = num2str(comps(i));
    r_file = ['results_meg/not_ld/avg_megsub0_rvals_100_',comp,'_feature3.mat']; 
    rVals = load(r_file,'rVals').rVals;
    
    num_trials = 4;
    curr = zeros(num_trials,1);
    for row = 1:num_trials
        total = 0;
        for col = 1:size(rVals,2)
            val = rVals(row,col);
            total = total + abs(val);
        end
        total = total / size(rVals,2);
        curr(row, 1) = total;
    end
    
    % h = bar(curr);
    % run prepExport.m
    
    % testAccuracy(curr(:, model));
    disp(['Mean r (for num comps=',comp,') = ',num2str(mean(curr(:,1)))])
end

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
stds = [];
num_trials = 4;
seNum = sqrt(num_trials);
comps = [64, 16, 64];
comps_ld = [16,16,16];
for idx = 3:5
    comp = num2str(comps(idx-2));
    comp_ld= num2str(comps_ld(idx-2));
    rnd_r_file = ['results_meg/rnd/rndimmegsub0_rvals_100_155_feature', num2str(idx), '.mat'];
    rndR = load(rnd_r_file,'rVals').rVals;

    indiv_r_file = ['results_meg/indiv/immegsub5_rvals_100_155_feature', num2str(idx), '.mat'];
    indivR = load(indiv_r_file,'rVals').rVals;
    
    avg_r_file = ['results_meg/avg_mod/average_immegsub5_rvals_100_155_feature', num2str(idx), '.mat'];
    avgR = load(avg_r_file,'rVals').rVals;
    
    mcca_r_file = ['results_meg/not_ld/imavg_megsub0_rvals_100_',comp,'_feature', num2str(idx), '.mat'];
%     if idx == 5
%         mcca_r_file = ['results_mcca/cmb/cmb_rvals_32_feature', num2str(idx), '.mat'];
%     end
    mccaR = load(mcca_r_file,'rVals').rVals;

    mccald_r_file = ['results_meg/ld/imavg_megldsub0_rvals_100_',comp_ld,'_feature', num2str(idx), '.mat'];
    mccaRld = load(mccald_r_file,'rVals').rVals;
    

    data = {rndR;indivR;avgR;mccaR;mccaRld};
    curr = zeros(num_trials,5);

    for i = 1:5
        currData = data{i,:};
        for row = 1:num_trials
            total = 0;
            for col = 1:size(currData,2)
                val = currData(row,col);
                total = total + abs(val);
            end
            total = total / size(currData,2);
            curr(row, i) = total;
        end
    end
    
    means = [means; [mean(curr(:,1)) mean(curr(:,2)) mean(curr(:,3)) mean(curr(:,4)) mean(curr(:,5))]];
    stds = [stds; [std(curr(:,1))/seNum std(curr(:,2))/seNum std(curr(:,3))/seNum std(curr(:,4))/seNum mean(curr(:,5))/seNum]];
end
 

x = categorical(["F0" "Spectrogram" "Aperiodicity"]);
h = bar(x,means);
hold on

[ngroups, nbars] = size(means);
pos = nan(nbars, ngroups);
for i = 1:nbars
    pos(i,:) = h(i).XEndPoints;
end
errorbar(pos',means,stds,'k','linestyle','none');

title('Mean R Value per Stimulus Feature'), xlabel('Stimulus Feature'), ylabel('R')
run prepExport.m
legend('Random Model','Individual Model','Average Model','MCCA Model', 'MCCA Model (Small λ)', '', '','','') 

%% evaluating individual components

names = ["F0", "Spectrogram", "Aperiodicity"];
seNum = sqrt(20);
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
    stds = zeros(16,1);
    for j = 1:16
        %testAccuracy(curr(:, model));
        disp(['Mean r for component ', num2str(j), ' [feature ', num2str(stimIdx), '] = ',num2str(mean(curr(:,j)))])
        means(j,1) = mean(curr(:,j));
        stds(j,1) = std(curr(:,j))/seNum ;
    end
    
    h = bar(means);
    hold on

    [ngroups, nbars] = size(means);
    pos = nan(nbars, ngroups);
    for i = 1:nbars
        pos(i,:) = h(i).XEndPoints;
    end
    errorbar(pos',means,stds,'k','linestyle','none');
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

% Loading prediction and stimulus data
predFile = ['results_mcca/cmb/cmb_pred_32_feature3.mat'];
stimFile = 'stim_input_files/dataStim_32.mat';
stim = load(stimFile,'stim').stim;

stimIdx=3;
num_audios = 20;
highValsMCCA = [];
allValsMCCA = [];

for i = 1:num_audios
    mccaPred = load(predFile,'predAll').predAll{1,i};
    feature = stim.data{stimIdx,i};
    feature = feature(1:size(mccaPred,1),:);
    
    highPred = [];
    highFeat = [];
    theNum = quantile(feature,0.5);
    for j = 1:size(mccaPred,1)
        if feature(j,1) > theNum
            highPred = [highPred; mccaPred(j)];
            highFeat = [highFeat; feature(j)];
        end
    end

    allMccaCorr = corrcoef(feature,mccaPred);
    allValsMCCA = [allValsMCCA allMccaCorr(1,2)];

    mccaCorr = corrcoef(highFeat,highPred);
    highValsMCCA = [highValsMCCA mccaCorr(1,2)];
end

disp(['Mean corr for high val MCCA = ',num2str(mean(highValsMCCA))])
disp(['Mean corr for all MCCA = ',num2str(mean(allValsMCCA))])

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

stimIdx = 4;
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

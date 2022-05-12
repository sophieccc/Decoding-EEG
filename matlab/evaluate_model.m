clear variables
close all

addpath libs/cnsp_utils
addpath libs/cnsp_utils/cnd
addpath libs/mTRF-Toolbox_v2/mtrf
addpath libs/NoiseTools
addpath eeglab2021.1
eeglab

%%
stimIdx = 5;

%% Compare r values

indiv_r_file = ['indiv_sub_results/sub1_rvals_32_feature', num2str(stimIdx), '.mat'];
indivR = load(indiv_r_file,'rVals').rVals;

mcca128_r_file = ['mcca_32_results/sub0_rvals_32_128_feature', num2str(stimIdx), '.mat'];
mccaR128 = load(mcca128_r_file,'rVals').rVals;

mcca_r_file = ['mcca_32_results/sub0_rvals_32_feature', num2str(stimIdx), '.mat'];
mccaR = load(mcca_r_file,'rVals').rVals;

mcca64_r_file = ['mcca_64_results/sub0_rvals_64_feature', num2str(stimIdx), '.mat'];
mcca64R = load(mcca64_r_file,'rVals').rVals;

avg_r_file = ['avg_sub_results/avg_sub1_rvals_32_feature', num2str(stimIdx), '.mat'];
avgR = load(avg_r_file,'rVals').rVals;

data = {indivR;mccaR128;mccaR;mcca64R;avgR};
curr = zeros(20,5);

for i = 1:size(data,1)
    currData = data{i,:};
    for row = 1:20
        total = 0;
        for col = 1:size(currData,2)
            val = currData(row,col);
            total = total + val;
        end
        total = total / size(currData,2);
        curr(row, i) = total;
    end
end

h = bar(curr);
set(h, {'DisplayName'}, {'Individual','MCCA32', 'MCCA32 (128)','MCCA64','Average'}')
run prepExport.m
legend() 

for model = 1:5
    %testAccuracy(curr(:, model));
    disp(['Mean r = ',num2str(mean(curr(:,model)))])
end


%%

% Loading mcca model data
modelFilename = ['mcca_32_results/sub0_model_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','sub0_model_32_feature', num2str(stimIdx), '.mat']);
mccaModel = load(modelFilename,'modelAll').modelAll;

% Loading avg model data
avgModelFilename = ['avg_sub_results/avg_model_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','avg_model_32_feature', num2str(stimIdx), '.mat']);
avgModel = load(avgModelFilename,'modelAll').modelAll;


mccaAvgModel = plotTRF(mccaModel);
plotGFP(mccaAvgModel);

avgAvgModel = plotTRF(avgModel);
plotGFP(avgAvgModel);

%%

% Loading mcca prediction data
mccaFile = ['mcca_32_results/sub0_pred_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','sub0_pred_32_feature', num2str(stimIdx), '.mat']);
mccaPred = load(mccaFile,'predAll').predAll{1};

% Loading avg prediction data
avgFile = ['avg_sub_results/avg_pred_32_feature', num2str(stimIdx), '.mat'];
disp(['Loading model data: ','avg_pred_32_feature', num2str(stimIdx), '.mat']);
avgPred = load(avgFile,'predAll').predAll{1};
avgPred = avgPred(1:size(mccaPred,1),:);

% Loading stimulus data
stimFile = 'dataStim_32.mat';
disp('Loading model data: dataStim_32.mat');
stim = load(stimFile,'stim').stim;
feature = stim.data{stimIdx,1};
feature = feature(1:size(mccaPred,1),:);

mccaCorr = corrcoef(feature,mccaPred)
avgCorr = corrcoef(feature,avgPred)

%%

function testAccuracy(rAll)
    disp(size(rAll))
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
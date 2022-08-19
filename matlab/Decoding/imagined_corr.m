clear variables
result_size = 100;
squeezeTime = 1;
predFile = 'results_meg/ld/combavg_megldsub0_pred_100_64_feature4.mat';

overall_top = [];
overall_bottom = [];
overall_pred_top = [];
overall_pred_bottom = [];

for i = 1:2
    [pred, stimData, onsets] = loadData(i, predFile);

    % (skipping final word for ease.)
    numWords = size(onsets,1) - 1;

    % Getting stimulus word correlations.
    stimCorr = getCorrMat(stimData, numWords, onsets, squeezeTime);

    % empty half of diagonal due to duplicate info.
    stimCorr2 = emptyDiagonal(stimCorr, numWords);

    % get top and bottom 100 correlations.
    [topWords, bottomWords] = getXCorrelations(result_size, numWords, stimCorr2);
    overall_top = [overall_top ;topWords];
    overall_bottom = [overall_bottom ;bottomWords];

    predCorr = getCorrMat(pred, numWords, onsets, squeezeTime);
    [topPred, bottomPred] = getMatchingPredictions(result_size, topWords,bottomWords, predCorr);
    overall_pred_top = [overall_pred_top ;topPred];
    overall_pred_bottom = [overall_pred_bottom ;bottomPred];
end

[topWords, bottomWords] = pickFinalWords(result_size, overall_top, overall_bottom);
[topPred, bottomPred] = getFinalPreds(result_size, topWords, bottomWords, overall_pred_top, overall_pred_bottom);

%% Plotting results.

sqn = sqrt(result_size);

% Stimulus vs Imagery
% means = [[mean(topWords(:,1)) mean(bottomWords(:,1))];[mean(topPred(:,1)) mean(bottomPred(:,1))]];
% stds = [[std(topWords(:,1))/sqn std(bottomWords(:,1))/sqn];[std(topPred(:,1))/sqn std(bottomPred(:,1))/sqn]];
% x = categorical(["Stimulus" "Imagined"]);

%imagData = load("imag_corrdata.mat",'stuff').stuff;

% Run on imagined data, note values, then fill in here when running on
% Listened data to get full bar chart.
imagTop =  0.6581; %imagData(6,1); 
imagBottom = 0.4163; %imagData(6,2); 
imagSeTop = 0.0286;%imagData(6,3); 
imagSeBottom = 0.0268;%imagData(6,4); 

% Stimulus vs Listened vs Imagery
stimulusMeans = [mean(topWords(:,1)) mean(bottomWords(:,1))];
predMeans = [mean(topPred(:,1)) mean(bottomPred(:,1))];
means = [stimulusMeans;[imagTop imagBottom];predMeans];

stimulusSes = [std(topWords(:,1))/sqn std(bottomWords(:,1))/sqn];
predSes = [std(topPred(:,1))/sqn std(bottomPred(:,1))/sqn];
stds = [stimulusSes;[imagSeTop imagSeBottom];predSes];

predMeans
predSes

means
stds

plotBars(means,stds);

%% IGNORE. Testing stuff to see if it works

% Set stimCorr2, stimCorr3, and predCorr values manually.
[topWords, bottomWords] = getXCorrelations(3, 6, stimCorr2);
[topWords2, bottomWords2] = getXCorrelations(3, 6, stimCorr3);
overall_top = [topWords ;topWords2];
overall_bottom = [bottomWords ;bottomWords2];

[topPred, bottomPred] = getMatchingPredictions(3, topWords, bottomWords, predCorr);
[topPred2, bottomPred2] = getMatchingPredictions(3, topWords2, bottomWords2, predCorr2);
overall_pred_top = [topPred ;topPred2];
overall_pred_bottom = [bottomPred ;bottomPred2];

[finalTopWords, finalBottomWords] = pickFinalWords(3, overall_top, overall_bottom);
[finalTopPred, finalBottomPred] = getFinalPreds(3, finalTopWords, finalBottomWords, overall_pred_top, overall_pred_bottom);

%% Functions.

function [corrMat] = getCorrMat(theData, numWords, onsets, squeezeTime)
    corrMat  = zeros(numWords,numWords);
    for i = 1:numWords
        for j = 1:numWords
            if i == j
                corrMat(i,j) = 0;
            else
                firstWord = theData(onsets(i):onsets(i+1),:);
                secondWord = theData(onsets(j):onsets(j+1),:);
           
                if squeezeTime == 1
                    firstWord = mean(firstWord,1);
                    secondWord = mean(secondWord,1);
                else
                    wordSize = min(size(firstWord,1), size(secondWord,1));
                    firstWord = firstWord(1:wordSize,:);
                    secondWord = secondWord(1:wordSize,:);
                end

                %if i == 47 && j ==40 && (size(onsets,1) < 60)
                if i == 32 && j ==18 && (size(onsets,1) < 60)
%                     plotSpectro(firstWord');
%                     plotSpectro(secondWord');
                end

                r = corrcoef(firstWord,secondWord);
                corrMat(i,j) = r(1,2);
            end
        end
    end
end

function [topWords, bottomWords] = getXCorrelations(num_res, numWords, stimCorr2)
    topWords = zeros(num_res,3);
    bottomWords = ones(num_res,3);
    for i = 1:numWords
        for j = 1:numWords
            % To avoid the blank diagonal.
            if i > j
                [miny,x] = min(topWords(:,1));
                [maxy,y] = max(bottomWords(:,1));
                if stimCorr2(i,j) > miny
                    topWords(x,1) = stimCorr2(i,j);
                    topWords(x,2) = i;
                    topWords(x,3) = j;
                end
                if abs(stimCorr2(i,j)) < maxy
                    bottomWords(y,1) = abs(stimCorr2(i,j));
                    bottomWords(y,2) = i;
                    bottomWords(y,3) = j;
                end
    
            end
        end
    end
end

function [pred, stimData, onsets] = loadData(num, predFile)
    clear [bottomWords, topWords, bottomPred, topPred];
    pred = load(predFile,'predAll').predAll{1,num};
    
    % Loading in stimulus data.
    stimData = load("stim_input_files/meg_audio_files100_4.mat", "stim").stim.data{4,num};
    stimData = stimData(1:size(pred,1),:);
    onsets = [];
    if num == 1
        onsetData = load("../onsets.mat", "poem1_onsets").poem1_onsets;
        onsets = find(onsetData);
    else
        onsetData = load("../onsets.mat", "poem2_onsets").poem2_onsets;
        onsets = find(onsetData);
    end
end

function stimCorr2 = emptyDiagonal(stimCorr, numWords)
    % empty half of diagonal due to duplicate info.
    stimCorr2 = stimCorr;
    for i = 1:numWords
        for j = 1:numWords
            if i < j
                stimCorr2(i,j) = 0;
            end
        end
    end
end

function [topPred, bottomPred] = getMatchingPredictions(result_size, topWords, bottomWords, predCorr)
    topPred = zeros(result_size,1);
    bottomPred = zeros(result_size,1);

    for j = 1:result_size
        topPred(j) = abs(predCorr(topWords(j,2), topWords(j,3)));
        bottomPred(j) = abs(predCorr(bottomWords(j,2), bottomWords(j,3)));
    end
end

function [topWords, bottomWords] = pickFinalWords(result_size, top, bottom)
    topWords = zeros(result_size,2);
    bottomWords = ones(result_size,2);
    for i = 1:(result_size * 2)
        [miny,x] = min(topWords(:,1));
        [maxy,y] = max(bottomWords(:,1));
        if top(i,1) > miny
            topWords(x,:,:,:) = [top(i,1),i];
        end
        if abs(bottom(i,1)) < maxy
            bottomWords(y,:,:,:) = [bottom(i,1),i];
        end
    end
end

function [topPred, bottomPred] = getFinalPreds(result_size, topWords, bottomWords, overallTop, overallBottom)
    topPred = zeros(result_size,1);
    bottomPred = zeros(result_size,1);

    for j = 1:result_size
        topPred(j) = overallTop(topWords(j,2));
        bottomPred(j) = overallBottom(bottomWords(j,2));
    end
end

function plotSpectro(data)
    figure;
    %hold on
    title('Word X: Y')
    xlabel('Time (Samples)')
    ylabel('Frequency Bands')
    imagesc(data)
    %run prepExport.m
end

function plotBars(means,stds)
    figure;
    x = categorical(["Stimulus" "Imagined" "Listened"]);
    h = bar(x,means);
    hold on
    
    [ngroups, nbars] = size(means);
    pos = nan(nbars, ngroups);
    for i = 1:nbars
        pos(i,:) = h(i).XEndPoints;
    end
    errorbar(pos',means,stds,'k','linestyle','none');
    
    title('Comparing Spectrograms Most and Least Correlated Words'), ylabel('Correlation')
    run prepExport.m
    legend('Top Correlated Words','Bottom Correlated Words','','')
end


    
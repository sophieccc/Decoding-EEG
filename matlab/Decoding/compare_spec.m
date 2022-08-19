% Comparing correlation of stimulus word with matching predicted word VS
% stimulus word with random predicted word.

predFile = 'results_meg/not_ld/imavg_megsub0_pred_100_32_feature4.mat';
stimFile = 'stim_input_files/meg_audio_files100_4.mat';

stim = load(stimFile,'stim').stim.data{4,1};
pred = load(predFile,'predAll').predAll{1,1};
onsets = find(poem1_onsets);
actualCorr = zeros(size(onsets,1),1);
randCorr = zeros(size(onsets,1),1);

for i = 1:size(onsets,1)
    endIndex = size(pred,1);
    if i~=size(onsets,1)
        endIndex = onsets(i+1);
    end
    firstStim = stim(onsets(i):endIndex,:);

    firstPred = pred(onsets(i):endIndex,:);
    randNum = randi([1 size(onsets,1)],1,1);
    endIndex = size(pred,1);
    if randNum~=size(onsets,1)
        endIndex = onsets(randNum+1);
    end
    secondPred = pred(onsets(randNum):endIndex,:);

    wordSize = min(size(firstPred,1), size(secondPred,1));
    firstPred = firstPred(1:wordSize,:);
    secondPred = secondPred(1:wordSize,:);
    firstStim = firstStim(1:wordSize,:);

    r = corrcoef(firstStim,firstPred);
    actualCorr(i) = r(1,2);
    r = corrcoef(firstStim,secondPred);
    randCorr(i) = r(1,2);
end

mean(actualCorr)
mean(randCorr)


stimFile = "stim_input_files/meg_audio_files100_4.mat";
predFile = "results_meg/not_ld/avg_megsub0_pred_100_64_feature3.mat";
outputFileName = "results_meg/not_ld/scale_avg_megsub0_pred_100_64_feature3.mat";
numTrials = 4;

origStim = load(stimFile).stim.data(3,:);

ls = [];
us = [];
for idx = 1:numTrials
    ls = [ls prctile(origStim{1,idx},5)];
    us = [us prctile(origStim{1,idx},95)];
end
l = mean(ls);
u = mean(us);

predAll = load(predFile).predAll;
for i = 1:numTrials
    B = rescale(predAll{1,i},l,u);
    predAll{1,i} = B;
end
save(outputFileName,"predAll");

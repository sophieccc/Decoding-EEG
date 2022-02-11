function recreateStim()

origStruct = load("dataStim_32.mat");
origStim = origStruct.stim;
newStim = struct(origStruct.stim);
stimData = origStim.data;

for i = 3:5
    for j = 1:size(stimData,2)
        name = "avg_sub1_pred_32_feature" + i + ".mat";
        feat = load(name);
        featData = feat.predAll;
        x = featData{j};
        y = cast(x,'double');
        y = resample(y, 200, 32);
        stimData{i,j} = y;
    end
end

newStim.data = stimData;
stim = newStim;
save("newStim32Avg.mat","stim");
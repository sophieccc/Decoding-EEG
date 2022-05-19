function recreateStim(templateStim, outputName, fs)

origStruct = load(templateStim);
origStim = origStruct.stim;
newStim = struct(origStruct.stim);
stimData = origStim.data;

for i = 3:5
    for j = 1:size(stimData,2)
        name = "results_mcca_32/cmb_pred_32_feature" + i + ".mat";
        feat = load(name);
        featData = feat.predAll;
        x = featData{j};
        y = cast(x,'double');
        y = resample(y, 200, fs);
        stimData{i,j} = y;
    end
end

newStim.data = stimData;
stim = newStim;
save(outputName,"stim");
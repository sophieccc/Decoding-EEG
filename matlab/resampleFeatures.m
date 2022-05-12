function resampleFeatures(filenames)

[~,c] = size(filenames);
if (c > 1)
    stim = combineStimFiles(filenames);
    stimData = stim.data;
else
    fileData = load(filenames);
    stim = fileData.stim;
    stimData = stim.data;
end

for i = 5:6
    for j = 1:size(stimData,2)
        x = stimData{i,j};
        y = cast(x,'double');
        y = resample(y, 64, 200);
        stimData{i,j} = y;
    end
end

stim.fs = 64.0;
stim.data = stimData;
save("newStim64.mat","stim");
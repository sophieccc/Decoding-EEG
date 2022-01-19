function [stim] = resampleFeatures(filename)

fileData = load(filename);
stim = fileData.stim
stimData = stim.data;

for i = 3:size(stimData,1)
    for j = 1:size(stimData,2)
        x = stimData{i,j};
        y = resample(x, 128, 200);
        stimData{i,j} = y;
    end
end

stim.data = stimData;
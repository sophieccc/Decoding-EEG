function [stim] = combineStimFiles(filenames)

fileData = load(filenames(1));
stim = fileData.stim;
stimData = stim.data;

for i = 2:size(filenames,2)
    currFile = load(filenames(i));
    stimData = [stimData currFile.stim.data];
end

stim.data = stimData;
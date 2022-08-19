% This file is used to take the MCCA components and put them into
% the appropriate CND format to work with the rest of the code.

% Parameters to change
mccaFile = "yy_pre_128.mat";
eegFile = "combined_pre_subs.mat";
sampleCNDFile = "../../datasets/LalorNatSpeech/dataCND/pre_dataSub1.mat";
outputFile = "subData_pre_128.mat";
numTrials = 20;

mccaData = load(mccaFile).yy;
eegData = load(eegFile, "eeg").eeg.data(1,:);
subData = cell(1, numTrials);

ind = 1;
for obs = 1:numTrials
    [obsLen, chan] = size(eegData{1,obs});
    subData{1,obs} = mccaData(ind:ind+obsLen-1, :);
    ind = ind + obsLen;
end

origEeg = load(sampleCNDFile).eeg;
eeg = struct(origStruct.eeg);
eeg.data = subData;
save(outputFile,"eeg");

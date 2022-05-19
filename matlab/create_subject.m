yy = load("yy_pre64_128.mat");
mccaData = yy.yy;

eegFile = load("combined_64_subs.mat", "eeg");
eegData = eegFile.eeg.data(1,:);

subData = cell(1, 20);
ind = 1;

for obs = 1:20
    [obsLen, chan] = size(eegData{1,obs});
    subData{1,obs} = mccaData(ind:ind+obsLen-1, :);
    ind = ind + obsLen;
end

origStruct = load("../../datasets/LalorNatSpeech/dataCND/pre_dataSub1.mat");
origEeg = origStruct.eeg;
eeg = struct(origStruct.eeg);
eeg.data = subData;
save("subData_pre64_128.mat","eeg");

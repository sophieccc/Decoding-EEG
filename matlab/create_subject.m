yy = load("yy.mat");
mccaData = yy.yy;

eeg = load("combined_subs.mat", "eeg");
eegData = eeg.eeg.data(1,:);

subData = cell(1, 20);
ind = 1;

for obs = 1:20
    [obsLen, chan] = size(eegData{1,obs});
    ind+obsLen
    subData{1,obs} = mccaData(ind:ind+obsLen-1, :);
    ind = ind + obsLen;
end

eegStruct = {};
eegStruct.data = subData;
save("subData.mat","subData");

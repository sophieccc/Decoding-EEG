yy = load("yy2.mat");
mccaData = yy.yy;

eegFile = load("combined_subs.mat", "eeg");
eegData = eegFile.eeg.data(1,:);

origStruct = load("../../datasets/LalorNatSpeech/dataCND/dataSub1.mat");
origEeg = origStruct.eeg;

for comp = 1:16
    subData = cell(1, 20);
    ind = 1;
    
    for obs = 1:20
        [obsLen, chan] = size(eegData{1,obs});
        subData{1,obs} = mccaData(ind:ind+obsLen-1, comp);
        ind = ind + obsLen;
    end
    
    eeg = struct(origStruct.eeg);
    eeg.data = subData;
    name = "comps/subData_comp" + num2str(comp) + ".mat";
    save(name,"eeg");
end


% creates an eeg subject that is the average of the 19 subjects 
avgData = cell(1,20);

for i = 1:20
    eegPreFilename = '../datasets/LalorNatSpeech/dataCND/pre_dataSub1.mat';
    load(eegPreFilename,'eeg')
    theMatrix = eeg.data{1,i};
    
    for j = 2:19
        eegPreFilename = ['../datasets/LalorNatSpeech/dataCND/pre_dataSub',num2str(j),'.mat'];
        load(eegPreFilename,'eeg')
        thisSub = eeg.data{1,i};

        prevLen = size(theMatrix,1);
        currLen = size(thisSub,1);
        minLen = min(prevLen,currLen);
        theMatrix = theMatrix(1:minLen,:);
        thisSub = thisSub(1:minLen,:);

        theMatrix = theMatrix + thisSub;
    end

    theMatrix = theMatrix / 19;
    avgData{1,i} = theMatrix;

end


origStruct = load("../datasets/LalorNatSpeech/dataCND/pre_dataSub1.mat");
origEeg = origStruct.eeg;
eeg = struct(origStruct.eeg);
eeg.data = avgData;
save("avgSubData.mat","eeg");

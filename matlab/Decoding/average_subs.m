% creates an eeg subject that is the average of the 19 subjects 

num_subs = 5;
num_trials = 20;
subNames = [2383, 2697, 2816, 2818, 2820];
outputFileName = "average_immegSubData.mat";

eegPreFilename = '../pre_imdataSubR2383.mat';
load(eegPreFilename,'eeg')

avgData = cell(1,num_trials);
for i = 1:num_trials
    theMatrix = eeg.data{1,i};
    for j = 2:num_subs
        eegPreFilename = ['../pre_imdataSubR', num2str(subNames(j)), '.mat'];
        load(eegPreFilename,'eeg')

        thisSub = eeg.data{1,i};
        prevLen = size(theMatrix,1);
        currLen = size(thisSub,1);
        minLen = min(prevLen,currLen);

        theMatrix = theMatrix(1:minLen,:);
        thisSub = thisSub(1:minLen,:);
        theMatrix = theMatrix + thisSub;
    end

    theMatrix = theMatrix / num_subs;
    avgData{1,i} = theMatrix;

end


origStruct = load(eegPreFilename);
origEeg = origStruct.eeg;
eeg = struct(origStruct.eeg);
eeg.data = avgData;
save(outputFileName,"eeg");

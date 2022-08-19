% This file combines all of the subjects trials into one big matrix.
% The trials are cut to the length of the shortest trial as they
% need to be the same length when used in the MCCA code.

% Parameters to change
stimFilename = 'datasets/LalorNatSpeech/dataCND/dataStim.mat';
outputFile = "Decoding/mcca/imcombined_meg_subs.mat";
numSubs = 5;
numTrials = 20;
prefix = "pre_imdataSub";
subs = ["R2383", "R2697", "R2816", "R2818", "R2820"]; %  for MEG

% prefix = "datasets/LalorNatSpeech/dataCND/pre_dataSub"; % for EEG
% subs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"];

load(stimFilename,'stim')
stimFeature = stim;
stimFeature.data = stimFeature.data(1,:);

data = cell(numSubs,numTrials);
for sub = 1:numSubs
    filename = prefix + subs(sub) + ".mat";
    fileData = load(filename);
    eeg = fileData.eeg;
    eegData = eeg.data;
    for obs = 1:numTrials
        minLen = 2745; % trial/error, need to find the shortest trial len
        eegData{obs} = double(eegData{obs}(1:minLen,:));

        [dim1, dim2] = cellfun(@size,eegData(:,obs),'uni',false);
        mat = cell2mat(eegData(:,obs));
        mat = mat2cell(mat,dim1{1,1}, dim2{1,1});
        data(sub,obs) = mat;
    end
end
eeg = {};
eeg.data = data;
save(outputFile,"eeg","-v7.3");

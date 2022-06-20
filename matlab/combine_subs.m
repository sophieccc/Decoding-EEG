
stimFilename = 'datasets/LalorNatSpeech/dataCND/dataStim.mat';
load(stimFilename,'stim')
stimFeature = stim;
stimFeature.data = stimFeature.data(1,:);
%     minSizes = zeros(20);
num_subs = 5;

data = cell(num_subs,20);
prefix = "pre_dataSub";
subs = ["R2383", "R2697", "R2816", "R2818", "R2820"]; % just for MEG
for sub = 1:num_subs
    filename = prefix + subs(sub) + ".mat";
    fileData = load(filename);
    eeg = fileData.eeg;
    eegData = eeg.data;
    for obs = 1:20

%             envLen = size(stimFeature.data{obs},1);
%             eegLen = size(eegData{obs},1);
%             minLen = min(envLen,eegLen);
%             if(minSizes(obs)==0 || minLen < minSizes(obs))
%                 minSizes(obs) = minLen;
%             else
%                 minLen = minSizes(obs);
%             end
        minLen = 2746; %10603; %5302; %21205;
        eegData{obs} = double(eegData{obs}(1:minLen,:));

        [dim1, dim2] = cellfun(@size,eegData(:,obs),'uni',false);
        mat = cell2mat(eegData(:,obs));
        mat = mat2cell(mat,dim1{1,1}, dim2{1,1});
        data(sub,obs) = mat;
    end
end
eeg = {};
eeg.data = data;
save("combined_meg_subs.mat","eeg","-v7.3");

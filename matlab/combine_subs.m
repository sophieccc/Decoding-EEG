
stimFilename = 'datasets/LalorNatSpeech/dataCND/dataStim.mat';
load(stimFilename,'stim')
stimFeature = stim;
stimFeature.data = stimFeature.data(1,:);
%     minSizes = zeros(20);

data = cell(19,20);
prefix = "datasets/LalorNatSpeech/dataCND/pre_dataSub";
for sub = 1:19
    filename = prefix + sub + ".mat";
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
        minLen = 10603; %5302; %21205;
        eegData{obs} = double(eegData{obs}(1:minLen,:));

        [dim1, dim2] = cellfun(@size,eegData(:,obs),'uni',false);
        mat = cell2mat(eegData(:,obs));
        mat = mat2cell(mat,dim1{1,1}, dim2{1,1});
        data(sub,obs) = mat;
    end
end
eeg = {};
eeg.data = data;
save("combined_64_subs.mat","eeg","-v7.3");

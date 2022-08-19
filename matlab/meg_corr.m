
% Loading in a subject's trial data. 'data' is a 1x20 cell.
% Each cell contains a listening trial (1-10 = poem 1, 11-20 = poem 2).
data = load("pre_dataSubR2820.mat", "eeg").eeg.data;

%% 
numTrials = 20;
numChannels = 155;

res = zeros(numTrials-1,numChannels);
for chan = 1:numChannels
    for tr = 2:numTrials
        % Assuming minimal preprocessing, i.e. different no. of samples
        minLen = min([size(data{1,1},1), size(data{1,tr},1)]);

        % Correlation between trial 1 and trial i.
        res(tr-1, chan) = corr(data{1,1}(1:minLen,chan), data{1,tr}(1:minLen,elec));
    end
end

% Averaging over channels.
res = mean(res,2);

% Comparing trial 1 (of poem 1) to the other poem 1 trials.
mean(res(1:(numTrials/2)-1))

% Comparing trial 1 (of poem 1) to the poem 2 trials.
mean(res(10:numTrials-1))
    
% For ease of evaluation when creating separate models for each band.

numTrials = 20;
numChannels = 32;
trialLen = 5302;

newRVals = zeros(numTrials,numChannels);
for i = 1:numTrials
    for j = 1:numChannels
        newRVals(i,j) = theRs{j}(i);
    end
end

rVals = newRVals;
name = "band_rvals_32_feature4.mat";
save(name,"rVals");


newPredVals = cell(1,numTrials);
for i = 1:numTrials
    curr = zeros(trialLen,numChannels);
    for j = 1:numChannels
        temp = thePreds{j,:};
        for x = 1:numTrials
            curr(:,j) = temp{1,i};
        end
    end
    newPredVals{1,i} = curr;
end

predAll = newPredVals;
outputFile = "band_pred_32_feature4.mat";
save(outputFile,"predAll");

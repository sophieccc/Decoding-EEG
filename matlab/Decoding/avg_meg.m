% avg before model (average MEG data)
numSamples = 2745;
megFile = load("mcca/imcombined_meg_subs.mat");
numTrials = 4;
numSubjects = 5;
numChannels = 155;

eeg = megFile.eeg;
newData = eeg.data(:,1:numTrials);

for sub = 1:numSubjects
    for x = 1:numTrials
        for i = 1:numSamples
            for j = 1:numChannels
                total = 0;
                first_num = (numSubjects * x) - numTrials;
                second_num = numSubjects * x;
                for x2 = first_num:second_num
                    curr = eeg.data{sub,x2};
                    total = total + curr(i,j);
                end
                newData{sub,x}(i,j) = total/numSubjects;
            end
        end
    end
end

eeg.data = newData;
name = "mcca/imcombined_meg_avg.mat";
save(name,"eeg");

%%
% avg after model (average predictions)

prefix = "results_meg/megsub0_pred_100_64_feature";
numSamples = 2746;
numTrials = 20;
outputFileName = "results_meg/avg_meg_100_64_feature";

for f = 3:5
    numFilters = 32;
    if f == 3
        numFilters = 1;
    end

    megFile = load(prefix + f + ".mat");
    preds = megFile.predAll;
    poem_1 = zeros(numSamples,numFilters);
    poem_2 = zeros(numSamples,numFilters);

    for i = 1:numSamples
        for j = 1:numFilters
            one_total = 0;
            two_total = 0;
            for x = 1:(numTrials/2)
                one_total = one_total + preds{1,x}(i,j);
            end
            for x = ((numTrials/2)+1):numTrials
                two_total = two_total + preds{1,x}(i,j);
            end
            poem_1(i,j) = one_total/(numTrials/2);
            poem_2(i,j) = two_total/(numTrials/2);
        end
    end

    predAll = cell(1,2);
    predAll{1,1} = poem_1;
    predAll{1,2} = poem_2;    
    name = outputFileName + f + ".mat";
    save(name,"predAll");
end

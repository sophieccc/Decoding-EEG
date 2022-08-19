% Redundant, now done using python. 

function reduceFilterNumber(filenames, numFilters, featureIdx, fs)

    [~,c] = size(filenames);
    if (c > 1)
        stim = combineStimFiles(filenames);
        stimData = stim.data;
    else
        fileData = load(filenames);
        stim = fileData.stim;
        stimData = stim.data;
    end
    
    [~, numObs] = size(stimData);
    feature = stimData{featureIdx,1};
    [~, origFilters] = size(feature);
    maxFreq = floor(fs / 2);
    stepSize = maxFreq / origFilters;
    [lower,~,upper]= greenwud(numFilters, 1, maxFreq, 0);
    
    for obs = 1:numObs
        obsData = stimData{featureIdx,obs};
        [numSamples, origFilters] = size(obsData);
        newMat = zeros(numSamples,numFilters);
        for i = 1:numSamples
            sampleData = obsData(i,:);
            for j = 1:numFilters
                bottomFilt = floor(lower(j) / stepSize)+1;
                topFilt = floor(upper(j) / stepSize)+1;
                if topFilt > origFilters
                    topFilt = origFilters;
                end
                vals = sampleData(bottomFilt:topFilt);
                totalFreq = sum(vals);
                newMat(i,j) = totalFreq / ((topFilt - bottomFilt)+1);
            end
        end
        obsData = newMat;
        stimData{featureIdx,obs} = obsData;
    end
    stim.data = stimData;
    save("dataStim.mat","stim", "-v7.3");
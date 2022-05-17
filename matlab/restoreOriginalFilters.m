function restoreOriginalFilters(filenames, numFils, featureIdx, fs, outputName)

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
    [numSamples, origFilters] = size(feature);
    newMat = zeros(numSamples,numFils);
    maxFreq = floor(fs / 2);
    stepSize = maxFreq / numFils;
    [lower,~,upper]= greenwud(origFilters, 1, maxFreq, 0);
    
    for obs = 1:numObs
        obsData = stimData{featureIdx,obs};
        [numSamples, origFilters] = size(obsData);
        for i = 1:numSamples
            sampleData = obsData(i,:);
            for j = 1:origFilters
                bottomFilt = floor(lower(j) / stepSize)+1;
                topFilt = floor(upper(j) / stepSize)+1;
                if topFilt > numFils
                    topFilt = numFils;
                end
                val = sampleData(j);
                for x = bottomFilt:topFilt
                    newMat(i,x) = val;
                end
            end
        end
        obsData = newMat;
        stimData{featureIdx,obs} = obsData;
    end
    stim.data = stimData;
    save(outputName,"stim", "-v7.3");
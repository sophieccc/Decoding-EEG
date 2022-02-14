function average_predictions(stimIdx, stimDim)

    fileEnd = "_pred_32_feature" + stimIdx + ".mat";
    predFilenames = ["sub1" + fileEnd, "sub2" + fileEnd, "sub3" + fileEnd, "sub4" + fileEnd, "sub5" + fileEnd];
    avgPred = cell(1,20);
    for i = 1:20
        vals = [];
        
        for curr = 1:size(predFilenames,2)
            currFile = predFilenames(curr);
            data = load(currFile);
            currPred = data.predAll{1,i};
            vals = cat(stimDim+1,vals,currPred);
        end
    
        currAvg = mean(vals,stimDim+1);
        avgPred{1,i} = currAvg;
    end
    
    name = "avg" + fileEnd;
    predAll = avgPred;
    save(name, 'predAll');
    
    end
    
newRVals = zeros(20,32);
for i = 1:20
    for j = 1:32
        newRVals(i,j) = theRs{j}(i);
    end
end

rVals =newRVals;
name = "band_rvals_32_feature4.mat";
save(name,"rVals");


newPredVals = cell(1,20);
for i = 1:20
    curr = zeros(5302,32);
    for j = 1:32
        hi = thePreds{j,:};
        for x = 1:20
            curr(:,j) = hi{1,i};
        end
    end
    newPredVals{1,i} = curr;
end

predAll = newPredVals;
name = "band_pred_32_feature4.mat";
save(name,"predAll");

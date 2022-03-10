
% How to run it

clear variables
numSubs = 19;
numObs = 20;
eegMat = constructEegMat(numSubs, numObs);
%eegMat = getSavedEegMat();

NPCS = 32;
fprintf('Running mCCA \b\n')
ccSub = 1;

numSamples = size(eegMat,1);
for sub=1:numSubs
    xAll = eegMat(ccSub, :, :);
    xAll = permute(xAll,[2,3,1]);
    fprintf('xAll:')
    size(xAll)
    [xpc,topcs{ccSub}] = alainPCA(xAll,NPCS);
    
    size(xpc)
    xx(ccSub, :, :, :) = [xpc];
    clear xAll
    
    ccSub = ccSub + 1;
end

[nSubj,nSamples,NPCS,nTrials] = size(xx);
xx = permute(xx,[2,4,3,1]);
xx = reshape(xx,size(xx,1)*size(xx,2),size(xx,3),size(xx,4));
C = nt_cov(xx(:,:)); 
A = nt_mcca(C,NPCS);
yy = xx(:,:)*A; % samples x mccs


yy = yy(:,1:NPCS); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 19 * 20 cells (inside each cell is time * channel)
function eegMat = constructEegMat(numSubs, numObs) 
    eeg = load("combined_subs.mat", "eeg");
    data = eeg.eeg.data;
    eegMat = [];
    
    for sub = 1:numSubs
        obsArray = [];
        for obs = 1:numObs
            %data(sub,:) = cellfun(@cell2mat, data(sub,:), 'Uniform', 0);
            hi = data{sub,obs};
            obsArray = cat(1, obsArray,hi);
        end
        eegMat = cat(3, eegMat, obsArray);
    end
    eegMat = permute(eegMat,[3,1,2]);
    eegStruct = {};
    eegStruct.data = eegMat;
    save("eegMat.mat","eegStruct","-v7.3");
end

function eegMat = getSavedEegMat() 
    eegStruct = load("eegMat.mat", "eegStruct");
    eegMat = eegStruct.eegStruct.data;
end

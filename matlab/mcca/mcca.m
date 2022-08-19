% Runs mcca on matrix of combined EEG subjects to get components.

clear variables

% Parameters to change
numSubs = 5;
numObs = 4;
NPCS = 128;
eegFilePath = "imcombined_meg_avg.mat";
mode = 1;

if mode == 1
    eegMat = constructEegMat(numSubs, numObs);
else
    eegMat = getSavedEegMat();
end
fprintf('Running mCCA \b\n')

for sub=1:numSubs
    xAll = eegMat(sub, :, :);
    xAll = permute(xAll,[2,3,1]);
    fprintf('xAll:')
    size(xAll)
    [xpc,topcs{sub}] = alainPCA(xAll,NPCS);
    
    size(xpc)
    xx(sub, :, :, :) = [xpc];
    clear xAll
end

[nSubj,nSamples,NPCS,nTrials] = size(xx);
xx = permute(xx,[2,4,3,1]);
xx = reshape(xx,size(xx,1)*size(xx,2),size(xx,3),size(xx,4));
disp('shape of xx = ');
disp(size(xx));

C = nt_cov(xx(:,:)); 
disp('shape of C = ');
disp(size(C));

A = nt_mcca(C,NPCS);
disp('shape of A = ');
disp(size(A));

yy = xx(:,:)*A; % samples x mccs
disp('shape of yy = ');
disp(size(yy));

yy = yy(:,1:NPCS); 
disp('shape of yy = ');
disp(size(yy));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numSubs * numObs cells (inside each cell is time * channel)
function eegMat = constructEegMat(numSubs, numObs, eegFilePath) 
    eeg = load(eegFilePath, "eeg");
    data = eeg.eeg.data;
    eegMat = [];
    
    for sub = 1:numSubs
        obsArray = [];
        for obs = 1:numObs
            temp = data{sub,obs};
            obsArray = cat(1, obsArray,temp);
        end
        eegMat = cat(3, eegMat, obsArray);
    end
    eegMat = permute(eegMat,[3,1,2]);
    eegStruct = {};
    eegStruct.data = eegMat;
    save("outputMat.mat","eegStruct","-v7.3");
end

function eegMat = getSavedEegMat() 
    eegStruct = load("outputMat.mat", "eegStruct");
    eegMat = eegStruct.eegStruct.data;
end

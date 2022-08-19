% Goes from 2 trials for each audio -> 1 trial for each audio.

directory = "results_meg/ld/";
predFile = 'avg_megldsub0_pred_100_128_feature4.mat';

pred = load(directory + predFile,'predAll').predAll;
avgPreds = cell(1,2);
avgPreds{1,1} = (pred{1,1} + pred{1,2})/2;
avgPreds{1,2} = (pred{1,3} + pred{1,4})/2;

predAll = avgPreds;
save(directory + "comb" + predFile, 'predAll');

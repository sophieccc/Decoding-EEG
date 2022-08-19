origStruct = load("stim_input_files/meg_audio_filesX.mat");
origStim = origStruct.stim;
newStim = struct(origStruct.stim);
stimData = origStim.data;
currFs = 100;
prefix = "results_meg/";
names = ["megsub5_pred_100_32_feature3.mat"; 
    "megsub5_pred_100_32_feature4.mat"; 
    "megsub5_pred_100_32_feature5.mat"
    ];
for i = 3:5
    for j = 1:size(stimData,2)
        name = prefix + names(i-2);
        feat = load(name);
        featData = feat.predAll;
        x = featData{j};
        % for getting best piece of audio
        %x = x(4400:4499,:);
        y = cast(x,'double');
        y = resample(y, 200, currFs);
        stimData{i,j} = y;
    end
end

newStim.data = stimData;
stim = newStim;
save("stim_output_files/megStim_indiv.mat","stim");
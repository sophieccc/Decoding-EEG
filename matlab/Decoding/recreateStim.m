% Combining the 3 predicted features into a new 'stim' to input into
% the vocoder.

megStimFile = "stim_input_files/meg_audio_files100_4.mat";
outputFileName = "stim_output_files/meg/scale_immcca.mat";
currFs = 100;
newFs = 200;
num_audios = 4;

prefix = "results_meg/not_ld/";
names = [
    "scale_imavg_megsub0_pred_100_32_feature3.mat"; 
    "scale_imavg_megsub0_pred_100_32_feature4.mat"; 
    "scale_imavg_megsub0_pred_100_32_feature5.mat"
    ];

origStim = load(megStimFile).stim;
newStim = struct(origStruct.stim);
stimData = origStim.data;

for i = 3:5
    for j = 1:num_audios
        name = prefix + names(i-2);
        feat = load(name);
        featData = feat.predAll;
        x = featData{j};
        y = cast(x,'double');
%         if i ==4
%             y(y<0) = 0;
%             y = y - (min(y));
%         end
        y = y(1:2745,:);
        y = resample(y, newFs, currFs);
        stimData{i,j} = y;
        stimData = stimData(:,1:num_audios);
    end
end

newStim.data = stimData;
newStim.fs = newFs;
stim = newStim;
save(outputFileName,"stim");
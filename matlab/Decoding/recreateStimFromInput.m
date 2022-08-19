
origStruct = load("stim_input_files/meg_audio_filesX.mat");
origStim = origStruct.stim;
newStim = struct(origStruct.stim);
stimData = origStim.data;

for i = 3:5
    for j = 1:size(stimData,2)
        x = stimData{i,j};
        y = cast(x,'double');
        if i ==3
            y = resample(y, 200, 100,0);
        else
            y = resample(y, 200, 100);
        end
        stimData{i,j} = y;
    end
end

newStim.data = stimData;
stim = newStim;
save("stim_output_files/meg_audio_filesX","stim");
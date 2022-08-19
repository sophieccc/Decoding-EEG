
fileData = load("stim_input_files/fixedStim32.mat");
stim = fileData.stim;
stimData = stim.data;

otherFileData = load("stim_input_files/audio_files.mat");
otherStim = otherFileData.stim;
otherStimData = otherStim.data;

    for j = 1:20
        x = otherStimData{3,j};
        y = cast(x,'double');
        %hd = getLPFilt(200,newFs/2);
        %for ii = 1:length(y)
        %    y(ii) = filtfilthd(hd,y(ii));
        %end
        y = resample(y, 32, 200, 0);
        stimData{3,j} = y;
    end

stim.data = stimData;
name = "newStim32";
save(name,"stim");
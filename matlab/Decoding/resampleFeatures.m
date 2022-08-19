
filenames = ["stim_input_files/meg_audio_filesX.mat"];
num_audios = 2;
newFS = 100;
[~,c] = size(filenames);
if (c > 1)
    stim = combineStimFiles(filenames);
    stimData = stim.data;
else
    fileData = load(filenames);
    stim = fileData.stim;
    stimData = stim.data;
end

for i = 3:5
    for j = 1:num_audios
        x = stimData{i,j};
        y = cast(x,'double');
        if i == 3
            %hd = getLPFilt(200,newFS/2);
            %for ii = 1:length(y)
            %    y(ii) = filtfilthd(hd,y(ii));
            %end
            y = resample(y, newFS, 200, 0);
        else
            y = resample(y, newFS, 200);
        end
        stimData{i,j} = y;
    end
end

stim.fs = newFS;
stim.data = stimData;
name = "meg_audio_filesXX.mat";
save(name,"stim");
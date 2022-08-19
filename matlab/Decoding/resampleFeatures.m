
filenames = ["stim_input_files/meg_audio_filesX2.mat"];
num_audios = 4;
oldFS = 200;
newFS = 100;
outputFileName = "meg_audio_files100_4.mat";

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
%         if i == 3
            %hd = getLPFilt(200,newFS/2);
            %for ii = 1:length(y)
            %    y(ii) = filtfilthd(hd,y(ii));
            %end
%         else
            y = resample(y, newFS, oldFS,0);
%         end
        stimData{i,j} = y;
    end
end

stim.fs = newFS;
stim.data = stimData;
save(outputFileName,"stim");
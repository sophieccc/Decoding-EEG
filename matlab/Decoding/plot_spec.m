% Example of plotting loaded spectrogram.
curr = stim.data{4,1};
h = imagesc(curr');
title('Spectrogram for an Example Audio'), xlabel('Time'), ylabel('Frequency Bands')
run prepExport.m
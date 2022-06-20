origStim = load("stim_input_files/audio_files.mat").stim.data(3,:);
ls = [];
us = [];
for idx = 1:20
    ls = [ls prctile(origStim{1,idx},10)];
    us = [us prctile(origStim{1,idx},90)];
end
l = mean(ls);
u = mean(us);

stim = load("stim_output_files/nonorm_mcca32Stim.mat").stim;
A = stim.data;
for i = 1:20
    B = rescale(A{3,i},l,u);
    A{3,i} = B;
end
stim.data = A;
save("scale_mcca32Stim.mat","stim");
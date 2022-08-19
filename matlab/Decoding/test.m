[x,fsSpeech] = audioread('audio5.wav');
env = abs(hilbert( x(:,1) ));
envRes = resample(env,128,fsSpeech);

%%

origStruct = load("stim_input_files/dataStim_32.mat");
origStim = origStruct.stim;
stimData = origStim.data;

for i = 3:5
    predAll = stimData(i,:);
    name = "eegfeature" + num2str(i);
    save(name,"predAll");
end



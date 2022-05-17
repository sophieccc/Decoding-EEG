data = load('pred.mat');
env = data.pred.';
%env = sum(env{1},2);
%env = stim.data{1,1};
env = resample(env,44100,128);
env(env < 0) = 0;
env = env/100;

x = rand(size(env)); % carrier
x = x*2-1; % centre between -1 and 1

recSpeech = x .* env; % reconstructed speech

figure;plot(env,'k','LineWidth',3)
hold on;plot(recSpeech,'r')

filename = 'test.wav';
audiowrite(filename,recSpeech,44100);

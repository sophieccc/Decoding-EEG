[x,fs]=audioread('../../Decoding-EEG/python/result_wavs/theStandard_1.wav');

[x1,fs1]=audioread('../../Decoding-EEG/python/result_wavs/scale_mcca32Stim_0.wav');


llen = min(length(x),length(x1)); % forcing the two signals to the same duration

x = x(1:llen);

x1 = x1(1:llen);


env = abs(hilbert(x1));
%env2 = resample(env,128,44100);
%hi = stim.data{1,1}(1:21209);
%corrcoef(env2,hi)

x1 = 10*x1.*env; % modulating the sound by the original amplitude envelope (and multiplying by 10, just to make this louder)

audiowrite('../../Decoding-EEG/python/result_wavs/scale_mcca32Stim_0_cheat.wav',x1,fs);

%%
[x1,fs1]=audioread('Decoding-EEG/python/result_wavs/mcca32Stim_0.wav');

env = abs(hilbert(x1));
%env2 = resample(env,128,44100);
%env2 = env2(1:22729);
%hi = stim.data{1,1};
%corrcoef(env2,hi)
%env(env<prctile(env,50)) = 0;

x1 = 2*x1.*env;

%audiowrite('../../Decoding-EEG/python/result_wavs/mcca32Stim_0_half_cheat.wav',x1,fs);

%%
[x,fs]=audioread('Decoding-EEG/python/result_wavs/theStandard_0.wav');

env = abs(hilbert(x));

modulatedNoise = 10*(rand(size(env))-0.5).*env;

audiowrite('../../Decoding-EEG/python/result_wavs/modulatedNoise.wav',modulatedNoise,fs);

%%
clear
load("results_mcca/cmb/sub0_pred_32_feature1.mat")
thePred = predAll{1};
%thePred(thePred<prctile(thePred,20)) = 0;
thePred = resample(thePred,44100,32);

[x1,fs]=audioread('../../Decoding-EEG/python/result_wavs/mcca32Stim_0.wav');
x1 = x1(1:size(thePred,1));
x1 = 10*x1.*thePred; 

audiowrite('../../Decoding-EEG/python/result_wavs/mcca32Stim_0_env.wav',x1,fs);
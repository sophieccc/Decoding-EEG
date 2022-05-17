[x,fsSpeech] = audioread('audio5.wav');
env = abs(hilbert( x(:,1) ));
envRes = resample(env,128,fsSpeech);
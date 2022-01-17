import numpy as np
import pyworld as pw
import scipy.io as sio
import soundfile as sf

import librosa
import librosa.display
import matplotlib.pyplot as plt

def main():
    # Getting audio data and sampling rate.
    x, fs = sf.read('data/audio5.wav')
    print(fs)
    print(x.shape)

    # Using stereo data so just use one side/column.
    x_left = x[:, 0]

    # Using a python wrapper for pyworld so array needs to be C-compatible.
    x_left = x_left.copy(order='C')

    # plt.figure(figsize=(15,4))
    # librosa.display.waveplot(x_left,sr=fs, max_points=50000.0)
    # plt.show()


    # Getting vocoder parameters.
    f0, spectrogram, aperiodicity = pw.wav2world(x_left, fs)
    time = np.sort(f0)
    vuv = (aperiodicity[:, 0] < 0.5).astype(np.float32)[:, None]

    # plt.plot(time, f0, color="red")
    # plt.show()
    # plt.imshow(np.transpose(spectrogram), origin="lower", aspect="auto")
    # plt.colorbar();
    # plt.show()

    f0 = np.log10(np.transpose([f0])+1)
    aperiodicity = np.log10(-np.transpose(np.transpose(aperiodicity))+1)
    spectrogram[spectrogram > 1] = 1
    np.power(spectrogram, .1) #compressing

    print('Output shape (f0):', f0.shape)
    print('Output shape (vuv):', vuv.shape)
    print('Output shape (aperiodicity):', aperiodicity.shape)
    print('Output shape (spectrogram):', spectrogram.shape)

    # Convert to CND/Matlab format
    features = {"f0": f0, "vuv": vuv, "spectrogram": spectrogram, "aperiodicity": aperiodicity}
    info = {"fs": 200, "features": features}
    sio.savemat('out.mat', mdict={'eeg': info})

if __name__ == '__main__':
    main()
import numpy as np
import pyworld as pw
import scipy.io as sio
import soundfile as sf

def analyse_audio(filename, compressed_data):
    # Getting audio data and sampling rate.
    x, fs = sf.read(filename)
    print('Audio sampling rate:',fs)
    print('Audio shape:',x.shape)

    # Using stereo data so just use one side/column.
    x_left = x[:, 0]

    # Using a python wrapper for pyworld so array needs to be C-compatible.
    x_left = x_left.copy(order='C')

    # Getting vocoder parameters.
    f0, spectrogram, aperiodicity = pw.wav2world(x_left, fs)
    vuv = (aperiodicity[:, 0] < 0.5).astype(np.float32)[:, None]
    f0 = np.transpose([f0])

    # Compressing parameters.
    if compressed_data:
        f0 = np.log10(f0+1)
        aperiodicity = np.log10(-np.transpose(np.transpose(aperiodicity))+1)
        spectrogram[spectrogram > 1] = 1
        np.power(spectrogram, .1)

    return fs, f0, vuv, aperiodicity, spectrogram

def synthesise_audio(filename, fs, compressed_data):
    mat = sio.loadmat(filename)
    f0 = mat['eeg']['data'][0, 0]['f0'][0, 0][:,0].copy(order='C')
    aperiodicity = mat['eeg']['data'][0, 0]['aperiodicity'][0, 0].copy(order='C')
    spectrogram = mat['eeg']['data'][0, 0]['spectrogram'][0, 0].copy(order='C')

    if compressed_data:
        f0 = 10 ** f0
        aperiodicity = 10 ** aperiodicity
        np.power(spectrogram, 10)

    audio = pw.synthesize(f0, spectrogram, aperiodicity, fs)
    return audio

def main():
    compressed_data = 0
    fs, f0, vuv, aperiodicity, spectrogram = analyse_audio('audio5.wav', compressed_data)
    print('Output shape (f0):', f0.shape)
    print('Output shape (aperiodicity):', aperiodicity.shape)
    print('Output shape (spectrogram):', spectrogram.shape)
    print('Output shape (vuv):', vuv.shape)

    # Convert to CND/Matlab format
    features = {"f0": f0, "vuv": vuv, "spectrogram": spectrogram, "aperiodicity": aperiodicity}
    info = {"fs": 200, "data": features}
    sio.savemat('out.mat', mdict={'eeg': info})

    audio = synthesise_audio('out.mat', fs, compressed_data)
    print('Output shape (audio):', audio.shape)
    sf.write('synthesised_audio_decompress.wav', audio, fs)

if __name__ == '__main__':
    main()
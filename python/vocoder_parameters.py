import os
import sys
import numpy as np
import pyworld as pw
import scipy.io as sio
import soundfile as sf

# The fs for all of our audio inputs is 44100.
input_fs = 44100

# The fs of the values returned by WORLD is always 200.
world_fs = 200

def analyse_audio(filename, compressed_data):
    # Getting audio data and sampling rate.
    audio, fs = sf.read(filename)

    # Using stereo data so just use one side/column.
    audio_left = audio[:, 0]

    # Using a python wrapper for pyworld so array needs to be C-compatible.
    audio_left = audio_left.copy(order='C')

    # Getting vocoder parameters.
    f0, spectrogram, aperiodicity = pw.wav2world(audio_left, fs)
    vuv = (aperiodicity[:, 0] < 0.5).astype(np.float32)[:, None]
    f0 = np.transpose([f0])

    # Compressing parameters.
    if compressed_data:
        f0 = np.log10(f0+1)
        aperiodicity = np.log10(-np.transpose(np.transpose(aperiodicity))+1)
        spectrogram[spectrogram > 1] = 1
        np.power(spectrogram, .1)

    return f0, spectrogram, aperiodicity, vuv

def format_features(f0, spectrogram, aperiodicity, vuv):
    features = [np.array(f0), np.array(spectrogram), np.array(aperiodicity), np.array(vuv)]
    data = np.zeros((4,), dtype=object)
    data[0] = features[0].astype('double')
    data[1] = features[1].astype('double')
    data[2] = features[2].astype('double')
    data[3] = features[3].astype('double')
    return data

def analyse_all_audios(directory_path, compressed_data):
    directory = os.fsencode(directory_path)
    data = []
    num_files = len(os.listdir(directory))

    for index in range(0, num_files):
        file_path = os.path.join(directory, os.listdir(directory)[index])
        f0, spectrogram, aperiodicity, vuv = analyse_audio(file_path, compressed_data)
        features = format_features(f0, spectrogram, aperiodicity, vuv)
        data.append(features)

    data_stim = sio.loadmat('data/dataStim.mat')
    stim_idx = data_stim['stim']['stimIdxs'][0][0].astype('double')
    cond_idx = data_stim['stim']['condIdxs'][0][0].astype('double')
    cond_names = data_stim['stim']['condNames'][0][0]
    stim_fs = data_stim['stim']['fs'][0][0]
    info = {"stimIdxs": stim_idx, "condIdxs": cond_idx, "condNames": cond_names, "fs": stim_fs}
    info["world_fs"] = world_fs
    info["data"] = np.concatenate((data_stim['stim']['data'][0][0][:, 0:num_files], np.transpose(data)))
    names = np.append(data_stim['stim']['names'][0][0][0], ["f0", "spectrogram", "aperiodicity", "vuv"])
    info["names"] = [names]
    sio.savemat("data/" + directory_path + '.mat', mdict={'stim': info})

def synthesise_audio(mat, fs, compressed_data):
    f0 = mat[2][:, 0].copy(order='C')
    spectrogram = mat[3].copy(order='C')
    aperiodicity = mat[4].copy(order='C')

    if compressed_data:
        f0 = 10 ** f0
        aperiodicity = 10 ** aperiodicity
        np.power(spectrogram, 10)

    audio = pw.synthesize(f0, spectrogram, aperiodicity, fs)
    return audio

def synthesise_all_audios(directory_path, compressed_data, fs):
    directory = os.fsencode(directory_path)
    mat = sio.loadmat("data/" + directory_path + '.mat')
    data = mat['stim']['data'][0,0]

    for index in range(0, len(os.listdir(directory))):
        audio = synthesise_audio(data[:,index], fs, compressed_data)
        sf.write('synthesised_' + str(index) + '.wav', audio, fs)

def main():
    compressed_data = 0
    path = ''
    args = sys.argv[1:]
    if len(args) > 0:
        path = args[0]
        if len(args) > 1:
            compressed_data = int(args[1])

    analyse_all_audios(path, compressed_data)

    synthesise_all_audios(path, compressed_data, input_fs)

if __name__ == '__main__':
    main()

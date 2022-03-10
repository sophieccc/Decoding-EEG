import os
import re
import sys
import numpy as np
import pyworld as pw
import scipy.io as sio
import soundfile as sf

import feature_utility as util

# The fs for all of our audio inputs is 44100.
INPUT_FS = 44100

# The fs of the values returned by WORLD is always 200.
WORLD_FS = 200


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
        f0 = np.log10(f0 + 1)
        aperiodicity = np.log10(-np.transpose(np.transpose(aperiodicity)) + 1)
        spectrogram[spectrogram > 1] = 1
        np.power(spectrogram, .1)

    return f0, spectrogram, aperiodicity, vuv


def format_features(f0, spectrogram, aperiodicity, vuv):
    features = [
        np.array(f0),
        np.array(spectrogram),
        np.array(aperiodicity),
        np.array(vuv)
    ]
    data = np.zeros((4, ), dtype=object)
    data[0] = features[0].astype('double')
    data[1] = features[1].astype('double')
    data[2] = features[2].astype('double')
    data[3] = features[3].astype('double')
    return data


def order_files(directory):
    files = os.listdir(directory)
    file_nums = []
    filepaths = []
    for audio_file in files:
        filepaths.append(os.path.join(directory, audio_file))
        audio_file = audio_file.decode("utf-8")
        audio_file = int(re.sub("[^0-9]", "", audio_file))
        file_nums.append(audio_file)
    file_dict = {file_nums[i]: filepaths[i] for i in range(len(file_nums))}
    sorted_dict = sorted(file_dict.items())
    sorted_paths = [pair[1] for pair in sorted_dict]
    return sorted_paths


def analyse_all_audios(directory_path, compressed_data):
    directory = os.fsencode(directory_path)
    data = []
    num_files = len(os.listdir(directory))
    sorted_paths = order_files(directory)

    for index in range(0, num_files):
        print("Analyzing audio ", str(index))
        file_path = sorted_paths[index]
        f0, spectrogram, aperiodicity, vuv = analyse_audio(
            file_path, compressed_data)
        features = format_features(f0, spectrogram, aperiodicity, vuv)
        data.append(features)

    data_stim = sio.loadmat('data/dataStim.mat')
    stim_idx = data_stim['stim']['stimIdxs'][0][0].astype('double')
    cond_idx = data_stim['stim']['condIdxs'][0][0].astype('double')
    cond_names = data_stim['stim']['condNames'][0][0]
    stim_fs = data_stim['stim']['fs'][0][0]
    info = {
        "stimIdxs": stim_idx,
        "condIdxs": cond_idx,
        "condNames": cond_names,
        "fs": stim_fs
    }
    info["world_fs"] = WORLD_FS
    info["data"] = np.concatenate(
        (data_stim['stim']['data'][0][0][:, 0:num_files], np.transpose(data)))
    names = np.append(data_stim['stim']['names'][0][0][0],
                      ["f0", "spectrogram", "aperiodicity", "vuv"])
    info["names"] = [names]

    info["data"] = util.reduce_filter_num(info, 32, 3, INPUT_FS)
    info["data"] = util.reduce_filter_num(info, 32, 4, INPUT_FS)
    sio.savemat("data/" + directory_path + '.mat', mdict={'stim': info})


def synthesise_audio(mat, fs, compressed_data):
    f0 = mat[2][:, 0].copy(order='C')
    spectrogram = util.restore_original_filters(mat[3], 1025, INPUT_FS)
    aperiodicity = util.restore_original_filters(mat[4], 1025, INPUT_FS)
    spectrogram = spectrogram.copy(order='C')
    aperiodicity = aperiodicity.copy(order='C')

    if compressed_data:
        f0 = 10**f0
        aperiodicity = 10**aperiodicity
        np.power(spectrogram, 10)

    audio = pw.synthesize(f0, spectrogram, aperiodicity, fs)
    return audio


def synthesise_all_audios(filepath, compressed_data, fs):
    mat = sio.loadmat("data/" + filepath + '.mat')
    data = mat['stim']['data'][0, 0]

    for index in range(0, 3):
        audio = synthesise_audio(data[:, index], fs, compressed_data)
        sf.write('mcca_synthesised_' + str(index) + '.wav', audio, fs)


def main():
    compressed_data = 0
    path = 'newStim32MCCA'
    args = sys.argv[1:]
    if len(args) > 0:
        path = args[0]
        if len(args) > 1:
            compressed_data = int(args[1])

    #analyse_all_audios(path, compressed_data)
    synthesise_all_audios(path, compressed_data, INPUT_FS)


if __name__ == '__main__':
    main()

import argparse
import os
import re

import numpy as np
import pyworld as pw
import scipy.io as sio
import soundfile as sf

import feature_utility as util

# The fs for all of our audio inputs is 44100.
INPUT_FS = 44100

# The fs of the values returned by WORLD is always 200.
WORLD_FS = 200

parser = argparse.ArgumentParser()
parser.add_argument( "-m",
                    "--mode",
                    help="Create vocoder params (1) or Use vocoder params (2)?",
                    type=int)
parser.add_argument("-c",
                    "--compressed",
                    help="Compress feature data?",
                    type=bool)
parser.add_argument("-f",
                    "--filepath",
                    help="Filepath for audio file or matlab features")
parser.add_argument("-n",
                    "--num_audios",
                    help="Number of trials to synthesise",
                    type=int)


def analyse_audio(filename, compressed_data):
    # Getting audio data and sampling rate.
    audio, fs = sf.read(filename)

    # Using stereo data so just use one side/column.
    audio_left = audio[:, 0]

    # Using a python wrapper for pyworld as array needs to be C-compatible.
    audio_left = audio_left.copy(order='C')

    # Getting vocoder parameters.
    f0, spectrogram, aperiodicity = pw.wav2world(audio_left, fs)
    vuv = (aperiodicity[:, 0] < 0.5).astype(np.float32)[:, None]
    f0 = np.transpose([f0])

    # Compressing parameters.
    if compressed_data:
        print('compressing stim data...')
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


def get_files(directory):
    files = os.listdir(directory)
    file_nums = []
    filepaths = []
    for audio_file in files:
        filepaths.append(os.path.join(directory, audio_file))
        audio_file = audio_file.decode("utf-8")
        file_nums.append(audio_file)
    file_dict = {file_nums[i]: filepaths[i] for i in range(len(file_nums))}
    sorted_dict = sorted(file_dict.items())
    sorted_paths = [
        pair[1].decode().rsplit('/', 1)[-1] for pair in sorted_dict
    ]
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

    data_stim = sio.loadmat(
        '../../CNSP-workshop2021_code/datasets/LalorNatSpeech/dataCND/dataStim.mat'
    )
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
    sio.savemat(directory_path + '.mat', mdict={'stim': info})


def synthesise_audio(mat, fs, compressed_data):
    f0 = mat[2][:, 0].copy(order='C')
    spectrogram = util.restore_original_filters(mat[3], 1025, INPUT_FS)
    #spectrogram = mat[3];
    #aperiodicity = mat[4];
    aperiodicity = util.restore_original_filters(mat[4], 1025, INPUT_FS)
    spectrogram = spectrogram.copy(order='C')
    aperiodicity = aperiodicity.copy(order='C')

    if compressed_data:
        f0 = 10**f0
        aperiodicity = 10**aperiodicity
        np.power(spectrogram, 10)

    audio = pw.synthesize(f0, spectrogram, aperiodicity, fs)
    return audio


def synthesise_all_audios(filepath, compressed_data, num_audios, fs):
    prefix = ""
    #prefix = "../../CNSP-workshop2021_code/CNSP_tutorial/stim_output_files/meg/improvements/"
    suffix = ""
    #suffix = ".mat"

    mat = sio.loadmat(prefix + filepath + suffix)
    data = mat['stim']['data'][0, 0]

    for index in range(0, num_audios):
        audio = synthesise_audio(data[:, index], fs, compressed_data)
        sf.write('result_wavs/meg/' + filepath + '_' + str(index) + '.wav',
                 audio, fs)


def main():
    path = 'audio_files'
    compressed_data = 0
    num_audios = 4
    mode = 1

    args = parser.parse_args()
    if args.filepath:
        path = args.filepath
    if args.compressed:
        compressed_data = 1
    if args.num_audios:
        num_audios = args.num_audios
    if args.mode:
        mode = args.mode

    if mode == 1:
        analyse_all_audios(path, compressed_data)
    else:
        synthesise_all_audios(path, compressed_data, num_audios, INPUT_FS)

        # directory = os.fsencode("../../CNSP-workshop2021_code/CNSP_tutorial/stim_output_files/meg/")
        # num_files = len(os.listdir(directory))
        # sorted_paths = get_files(directory)

        # for index in range(0, num_files):
        #     print("Analyzing audio ", str(index))
        #     file_path = sorted_paths[index]
        #     synthesise_all_audios(file_path, compressed_data, INPUT_FS)


if __name__ == '__main__':
    main()

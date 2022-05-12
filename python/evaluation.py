import soundfile as sf
#from sti import stiFromAudio, readwav
from pystoi import stoi
import pysepm


def main():
    files = [
        "avg32_synthesised_2.wav", "indiv32_synthesised_2.wav",
        "mcca32_128_synthesised_2.wav", "mcca32_synthesised_2.wav",
        "mcca64_synthesised_2.wav"
    ]

    for f in files:
        clean, fs = sf.read('audio_files/audio3.wav')
        clean = clean[:, 0]

        denoised, fs = sf.read('result_wavs/' + f)
        denoised = denoised[0:len(clean)]
        clean = clean[0:len(denoised)]

        # Clean and den should have the same length, and be 1D
        res = stoi(clean, denoised, fs, extended=False)
        print('data for ' + f + ": ")
        print('stoi: ' + str(res))
        print('Frequency-weighted Segmental SNR: ' +
              str(pysepm.fwSNRseg(clean, denoised, fs)))
        print('Segmental Signal-to-Noise Ratio: ' +
              str(pysepm.SNRseg(clean, denoised, fs)))
        print('Log-likelihood Ratio : ' + str(pysepm.llr(clean, denoised, fs)))
        print('Cepstrum Distance Objective Speech Quality Measure : ' +
              str(pysepm.cepstrum_distance(clean, denoised, fs)))
        print('Coherence and speech intelligibility index  : ' +
              str(pysepm.csii(clean, denoised, fs)))

    # # read audio
    # refAudio, refRate = readwav('audio_files/audio1.wav')
    # #degrAudio, degrRate = readwav('fs32/synthesised_0.wav')
    # degrAudio, degrRate = readwav('audio_files/audio1.wav')

    # # calculate the STI. Visually verify console output.
    # stis = stiFromAudio(refAudio, degrAudio, refRate, name='eval1.wav')

    # print("Test Result:")

    # # test result
    # if abs(stis - 0.63) < 0.002:
    #     print("OK")
    #     return 0
    # else:
    #     print("FAILED")
    #     return 1


if __name__ == '__main__':
    main()

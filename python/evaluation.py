import soundfile as sf
import pysepm
import numpy


def main():
    #files = ["rnd_mcca32Stim_0.wav", "indiv32Stim_0.wav", "avg32Stim_0.wav", "nonorm_mcca32Stim_0.wav", "mcca32Stim_0.wav", "theStandard_0.wav"] 
    #files = [ "random",  "individual", "average", "mcca_spec16", "mcca_spec16_normf", "mcca_ld", "mcca_ld_normf"]
    #files = [ "immcca_realf0",  "immcca_realsp", "immcca_realap", "env_immcca", "two_immcca", "scale_immcca", "scalef0_immcca", "env_scale_immcca"]
    files = ["env_scale_immcca", "env_scale_mcca"]

    all_stoi = []
    all_fwSNR = []
    all_cd = []
    num_audios = 4;
    for f in files:
        stoi = 0
        fwSNR = 0
        cd = 0
        for i in range(num_audios):
            audio_num = 1
            if i >=2:
                audio_num = 2
            clean, fs = sf.read('meg_audio_files/poem-' + str(audio_num) + '.wav')
            #clean, fs = sf.read('audio_files/audio' + str(audio_num) + '.wav')
            clean = clean[:, 0]

            recon_num = i

            if f == "two_immcca":
                recon_num = audio_num - 1
            denoised, fs = sf.read('result_wavs/meg/' + f + '.mat_' +str(recon_num) + '.wav')
            #denoised, fs = sf.read('result_wavs/' + f)
            denoised = denoised[0:len(clean)]
            clean = clean[0:len(denoised)]

            # Clean and then should have the same length, and be 1D
            stoi += pysepm.stoi(clean, denoised, fs, extended=False)
            fwSNR += pysepm.fwSNRseg(clean, denoised, fs)
            cd += pysepm.cepstrum_distance(clean, denoised, fs)

        stoi = stoi/num_audios
        all_stoi.append(stoi)
        fwSNR = fwSNR/num_audios
        all_fwSNR.append(fwSNR)
        cd = cd/num_audios
        all_cd.append(cd)

        print('data for ' + f + ": ")
        print('stoi: ' + str(stoi))
        print('Frequency-weighted Segmental SNR: ' + str(fwSNR))
        print('Cepstrum Distance Objective Speech Quality Measure : ' + str(cd))

    print("[[" + " ".join(str(x) for x in all_stoi) + "]; [" + " ".join(str(x) for x in all_fwSNR) + "]; [" + " ".join(str(x) for x in all_cd) + "]]")
    A = numpy.array([all_stoi, all_fwSNR, all_cd])
    numpy.savetxt("imaginedimprovedmeg.csv", A, fmt='%.5f')


if __name__ == '__main__':
    main()

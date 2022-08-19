# Decoding-EEG

## Python

The main two files are vocoder.py and evaluation.py

vocoder.py: Add "-m 1' to retrieve vocoder parameters from audios and use '-m 2' to reconstruct audios using vocoder parameters.

evaluation.py: Calculates the average stoi, fwSNRseg, and CD across trials for each audio included.

## Matlab

Once you have a .mat file containing your stimulus features (vocoder parameters) in CND format and EEG data in CND format, you can do the following:

1. Preprocess the EEG data using vocoder_reconstruction.m
2. Resample the stimulus features using resampleFeatures.m
3. Combine the preprocessed EEG data into one matrix using combine_subs.m
4. Use mcca.m to get MCCA components for the EEG matrix
5. Use create_subject to create a new CND format subject out of these components
6. Create a model using this subject with vocoder_reconstruction.m
7. Analyse the results for this model using evaluate_model.m

Parameter values will need to be changed to match the filenames, number of subjects, channels, trials, fs, min trial length, etc. of your data.


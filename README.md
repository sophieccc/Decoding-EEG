# Decoding-EEG

## Python

The main two files are vocoder.py and evaluation.py

vocoder.py: Add "-m 1' to retrieve vocoder parameters from audios and use '-m 2' to reconstruct audios using vocoder parameters.

evaluation.py: Calculates the average stoi, fwSNRseg, and CD across trials for each audio included.

## Matlab

Once you have a .mat file containing your stimulus features (vocoder parameters) in CND format and brain data in CND format, you can do the following:

1. If not already done: for MEG data, there is an extra step using meg.m to get the data in CND format
2. If not already done: Preprocess the brain data using vocoder_reconstruction.m
3. Resample the stimulus features using resampleFeatures.m
4. Combine the preprocessed brain data into one matrix using combine_subs.m
5. Use mcca.m to get MCCA components for the brain data matrix
6. Use create_subject to create a new CND format subject out of these components
7. Create a model using this subject with vocoder_reconstruction.m
8. Analyse the results for this model using evaluate_model.m

Parameter values will need to be changed to match the filenames, number of subjects, channels, trials, fs, min trial length, etc. of your data.


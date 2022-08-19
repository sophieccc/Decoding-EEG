# Decoding-EEG

## Python

The main two files are vocoder.py and evaluation.py

vocoder.py: Add "-m 1' to retrieve vocoder parameters from audios and use '-m 2' to reconstruct audios using vocoder parameters.

evaluation.py: Calculates the average stoi, fwSNRseg, and CD across trials for each audio included.
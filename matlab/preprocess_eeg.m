clear all
close all

addpath libs/cnsp_utils
addpath libs/cnsp_utils/cnd
addpath libs/mTRF-Toolbox_v2/mtrf
addpath libs/NoiseTools
addpath eeglab2021.1
eeglab


%% Parameters - Natural speech listening experiment
reRefType = 'Avg'; % or 'Mastoids'
bandpassFilterRange = [1,8]; % Hz (indicate 0 to avoid running the low-pass
                          % or high-pass filters or both)
                          % e.g., [0,8] will apply only a low-pass filter
                          % at 8 Hz
% changed by sophie
downFs = 32; % Hz. *** fs/downFs must be an integer value ***

stimIdx = 4; % 1: env; 2: word onset; 3: f0; 4: sp; 5: ap; 6: vuv;

%% Preprocess EEG - Natural speech listening experiment
% Loading EEG data
eegFilename = 'mcca/subData.mat';
load(eegFilename,'eeg')
eeg = cndNewOp(eeg,'Load'); % Saving the processing pipeline in the eeg struct

% Filtering - LPF (low-pass filter)
if bandpassFilterRange(2) > 0
    hd = getLPFilt(eeg.fs,bandpassFilterRange(2));
    
    % Filtering each trial/run with a for loop
%         for ii = 1:length(eeg.data)
%             eeg.data{ii} = filtfilthd(hd,eeg.data{ii});
%         end
    
    % Filtering each trial/run with a cellfun statement
    eeg.data = cellfun(@(x) filtfilthd(hd,x),eeg.data,'UniformOutput',false);
    
    % Filtering external channels
    if isfield(eeg,'extChan')
        for extIdx = 1:length(eeg.extChan)
            eeg.extChan{extIdx}.data = cellfun(@(x) filtfilthd(hd,x),eeg.extChan{extIdx}.data,'UniformOutput',false);
        end
    end
    
    eeg = cndNewOp(eeg,'LPF');
end

% Downsampling EEG and external channels
 eeg = cndDownsample(eeg,downFs);
 % stim{1} = cndDownsample(stim{1},downFs);

% Filtering - HPF (high-pass filter)
if bandpassFilterRange(1) > 0 
    hd = getHPFilt(eeg.fs,bandpassFilterRange(1));
    
    % Filtering EEG data
    eeg.data = cellfun(@(x) filtfilthd(hd,x),eeg.data,'UniformOutput',false);
    
    % Filtering external channels
    if isfield(eeg,'extChan')
        for extIdx = 1:length(eeg.extChan)
            eeg.extChan{extIdx}.data = cellfun(@(x) filtfilthd(hd,x),eeg.extChan{extIdx}.data,'UniformOutput',false);
        end  
    end
    
    eeg = cndNewOp(eeg,'HPF');
end

% Re-referencing EEG data
eeg = cndReref(eeg,reRefType);

% Removing initial padding (specific to this dataset)
if isfield(eeg,'paddingStartSample')
    for tr = 1:length(eeg.data)
        eeg.data{tr} = eeg.data{tr}(eeg.paddingStartSample,:);
        for extIdx = 1:length(eeg.extChan)
            eeg.extChan{extIdx}.data = eeg.extChan{extIdx}.data{tr}(1+eeg.paddingStartSample,:);
        end
    end
end

% Saving preprocessed data
eegPreFilename = 'pre_subData.mat';
save(eegPreFilename,'eeg')
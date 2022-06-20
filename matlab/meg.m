%%

data_file = sprintf('MEG-PoeMusic-Imagery/MEGraw/R2820/R2820_sns.sqd');
data = ft_read_data(data_file);
data = data';

%% Take MEG data and separate it into trials properly.

modes = ["listening"]; % "imagery"];
files = ["R2820", "R2818", "R2816", "R2697", "R2383"];
channels = [183, 183, 183, 183, 173];

for m = 1:size(modes,2)
    mode = modes(m);
    for i = 1:size(files,2)
        file = files(i);
        order_file = strcat("MEG-PoeMusic-Imagery/MEG-PoeMusic-Imagery/Log/",file,"/",file,".csv");
        order = get_order(order_file, mode);

        data_file = sprintf('MEG-PoeMusic-Imagery/MEGraw/%s/%s_sns.sqd', file, file);
        orig_data = ft_read_data(data_file);
        orig_data = orig_data';
        
        times = get_times(orig_data, channels(i));
        the_data = orig_data(:,1:157);
        data = split_data(the_data, times, order);

        data = data';
        save_file = strcat("dataSub", file, ".mat");
        save(save_file,'data');
    end
    
end

%% Take the data and do some MEG-specific preprocessing.
files = ["R2820", "R2818", "R2816", "R2697", "R2383"];
load("datasets/LalorNatSpeech/dataCND/dataSub1",'eeg')
old_eeg = eeg;
for i = 1:size(files,2)
    filename = strcat("dataSub", files(i), ".mat");
    load(filename,'data')
    data = remove_channels(data);
    convert_to_eeg_form(data, old_eeg, files(i));
end

%% Dealing with stimulus files (2 trials -> 20 trials). 

load("CNSP_tutorial/stim_input_files/meg_audio_files.mat",'stim')
oldData = stim.data;
newData = cell(6,20);

for i = 1:10
    for j = 1:6
        newData{j,i} = oldData{j,1};
    end
end
for i = 11:20
    for j = 1:6
        newData{j,i} = oldData{j,2};
    end
end

stim.data = newData;
stim.fs = 200;

save("CNSP_tutorial/stim_input_files/meg_audio_filesX.mat",'stim');
%%
function [data] = remove_channels(data)
    for i = 1:size(data,2)
        curr_cell = data{1,i};
        curr_cell(:,86) = [];
        curr_cell(:,56) = [];
        data{1,i} = curr_cell;
    end
end

function convert_to_eeg_form(data, old_eeg, file)
    eeg = struct(old_eeg);
    eeg.dataType = "MEG";
    eeg.fs = 1000;
    eeg.data = data;
    eeg = rmfield(eeg,"deviceName");
    eeg = rmfield(eeg,"extChan");
    eeg = rmfield(eeg,"chanlocs");
    filename = strcat("dataSub", file, ".mat");
    save(filename,'eeg');
end

function [order] = get_order(file, mode)
    t = readtable(file);
    poem_1 = zeros(10,1);
    poem_2 = zeros(10,1);
    p1_ind = 1;
    p2_ind = 1;
    for i = 1:size(t,1)
        if t{i,3} == mode & t{i,2} == "poem-1"
            poem_1(p1_ind) = i;
            p1_ind = p1_ind + 1;
        elseif t{i,3} == mode & t{i,2} == "poem-2"
            poem_2(p2_ind) = i;
            p2_ind = p2_ind + 1;
        end
    end
    order = [poem_1 poem_2];
    order=order(:);
end

function [resps] = split_data(data, times, order)
    resps = cell(20,1);
    for i = 1:size(order,1)
        first = times(order(i));
        second = size(data,1);
        if order(i) ~= size(times,1)
            second = times(order(i)+1);
        end
        curr_data = data(first:second,:);
        resps{i} = curr_data;
    end
end

function [times] = get_times(data, channel)
    x = data(:,channel);
    x = abs(x-median(x));
    x(x~=0) = 1;
    x = diff(x);
    x(x<0) = 0;
    times = find(x);
end




% This file is used to take each component from the first n MCCA components
% and create a new 'subject' using just that component.

% Parameters to change
n = 16;
num_trials = 20;
subDataFile = "subData_pre_128.mat";

subData = load(subDataFile).eeg.data;
for comp = 1:n
    compData = cell(1, num_trials);
    
    for obs = 1:num_trials
        compData{1,obs} = subData{1, obs}(:,comp);
    end
    
    eeg = struct(subDataFile.eeg);
    eeg.data = compData;
    name = "comps/subData_comp" + num2str(comp) + ".mat";
    save(name,"eeg");
end


subDataFile = load("subData_pre_128.mat");
subData = subDataFile.eeg.data;

for comp = 1:16
    compData = cell(1, 20);
    
    for obs = 1:20
        compData{1,obs} = subData{1, obs}(:,comp);
    end
    
    eeg = struct(subDataFile.eeg);
    eeg.data = compData;
    name = "comps/subData_comp" + num2str(comp) + ".mat";
    save(name,"eeg");
end


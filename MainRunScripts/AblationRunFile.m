%% Ablation Resting Recordings Main Runfile
% Input subject IDs as they correspond to collected data. This code will
% call run function without stimulation argin and process/output all
% relevant processed data. Processes the Resting data after removing each
% region from the model one at a time
allpts = {};
task = 'Resting';


for pat = 1:length(allpts)
    pt = allpts{pat};
    runfunctionablation(pt,task);
end

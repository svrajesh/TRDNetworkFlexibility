%% Resting State Recordings Main Runfile
% Input subject IDs as they correspond to collected data. This code will
% call run function without stimulation argin and process/output all
% relevant processed data.
allpts = {};
task = 'Resting';
for pat = 1:length(allpts)
    pt = allpts{pat};
    runfunhowction(pt,task);
end

%% Step03 -- Main Analysis Step
% This function runs the pMEM and ELA. Set the variable pt before running
% redo_PMEM and redo_Energy are toggleable in case only some parts need to
% be rerun, but generally safer to set both to 1 when doing full scale
% analysis

% parallel enabled (in nested level)

% Region Labels as encoded in region file
%   1:    Left vmPFC -- Gyrus Rectus, mOFC
%   2:    Left dmPFC -- Primarily SMA/SFG
%   3:    Left vlPFC -- lOFC/Pars/IFG
%   4:    Left dlPFC -- most of MFG, some SFG/IFG
%   5:    Left dACC
%   6:    Left MTL (Hipp, Amy)

%   7:    Right vmPFC -- Gyrus Rectus, mOFC
%   8:    Right dmPFC -- Primarily SMA/SFG
%   9:    Right vlPFC -- lOFC/Pars/IFG
%  10:    Right dlPFC -- most of MFG, some SFG/IFG
%  11:    Right dACC
%  12:    Right MTL (Hipp, Amy)

%   0: WM, Insula, Other

% method will be 'weighted'

%% PMEM analysis
% redo_PMEM = 1; %% CHANGE IF NOT RERUNNING
p = gcp("nocreate");
if isempty(p)
    parpool(32);
end
% if redo_PMEM
rootpath = betterLoadSave('USER_DIR');
switch method
      case {'weighted'}
        load(betterLoadSave(['USER_DIR/ChannelSelection/',task,'/',pt,'/final_recon.mat']));
end
% Current regions:
% exclude: MTL; dmPFC
allreg = 1:12; %1:14
keepreg = [1,3,4,5,7,9,10,11]; %[1,2,4,5,6,8,9,11,12,13];
[~,keepreg] = intersect(allreg,keepreg);
for ablatedregion = 1:length(keepreg)
    truekeepreg = keepreg;
    truekeepreg(ablatedregion) = [];
    final_selections_use = final_selections;
    for b = 1:length(final_selections)
        final_selections_use(b).channels = final_selections(b).channels(truekeepreg,:);
    end
    
    % % Create directories if they don't exist
    if ~isfolder(rootpath)
        mkdir(rootpath);
    end
    

    thresholds = 50;% only done at median threshold for this analysis
    for t = 1:length(thresholds)
        pMEM_model_fits(rootpath,pt,task,thresholds(t),final_selections_use,method,ablatedregion);
    end
    % end
end
%% Energy Landscape Analysis

% redo_Energy = 1;
% if redo_Energy
p = gcp("nocreate");
if isempty(p)
    parpool(32);
end


pmemfiles = dir(betterLoadSave(['USER_DIR/Analysis/Ablation/',task,'/',pt,'/pMEM/*/*.mat']));
pmemfiles = pmemfiles(~startsWith({pmemfiles.name},'.'));
pmemfiles = pmemfiles(~endsWith({pmemfiles.name},'Accuracies.mat'));
Nfiles = length(pmemfiles);
for fl = startIndex:Nfiles
    %tic
    flname = fullfile(pmemfiles(fl).folder,pmemfiles(fl).name);
    analyze_energy_landscape(flname);

end





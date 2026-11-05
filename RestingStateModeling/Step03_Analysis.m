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

load(betterLoadSave(['USER_DIR/ChannelSelection/',task,'/',pt,'/final_recon.mat']));

% Current regions:
% exclude: MTL; dmPFC
allreg = 1:12; %1:14
keepreg = [1,3,4,5,7,9,10,11];
[~,keepreg] = intersect(allreg,keepreg);
for b = 1:length(final_selections)
    final_selections(b).channels = final_selections(b).channels(keepreg,:);
end

% % Create directories if they don't exist
if ~isfolder(rootpath)
    mkdir(rootpath);
end

% binarization windows: given array is in seconds
thresholds = [25,30,35,37.5,40,45,50,55,60,62.5,65,70,75];
for t = 1:length(thresholds)
    pMEM_model_fits(rootpath,pt,task,thresholds(t),final_selections,'weighted',0);
end

%% Energy Landscape Analysis
% built in a logging feature in case this dies, as it was early in the
% project, but this is a waste now.


p = gcp("nocreate");
if isempty(p)
    parpool(32);
end


pmemfiles = dir(betterLoadSave(['USER_DIR/Analysis/',task,'/',pt,'/pMEM/*/*.mat']));
pmemfiles = pmemfiles(~startsWith({pmemfiles.name},'.'));
pmemfiles = pmemfiles(~endsWith({pmemfiles.name},'Accuracies.mat'));
Nfiles = length(pmemfiles);
for fl = startIndex:Nfiles
    flname = fullfile(pmemfiles(fl).folder,pmemfiles(fl).name);
    analyze_energy_landscape(flname);

end





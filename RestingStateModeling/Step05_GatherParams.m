%% Step05 -- Gathering params between energy and param files and adding a few things
% This technically could be wrappered in previous files but will keep this
% structure for clarity for now

%% Init folders
pMEMfolder = dir(betterLoadSave([USER_DIR,'/Analysis/',task,'/',pt,'/pMEM/*/*_result.mat']));
energyfolder = dir(betterLoadSave([USER_DIR,'/Analysis/',task,'/',pt,'/Energy/*/*/*.mat']));
pMEMfolder = pMEMfolder(~startsWith({pMEMfolder.name},'.'));
energyfolder = energyfolder(~startsWith({energyfolder.name},'.'));

paramfolder = betterLoadSave([USER_DIR,'/Analysis/',task,'/',pt,'/Params/']);

%% Loop folders
for i = 1:length(pMEMfolder)
    filename = pMEMfolder(i).name;
    foldername = pMEMfolder(i).folder;
    fullpmemfilename = fullfile(foldername,filename);
    bandname = extractBefore(filename,'_result');
    window = extractAfter(foldername,'features_');
    energyfilelist = energyfolder(strcmp(extractAfter({energyfolder.folder},'features_'),[window,'/',bandname]));
    
    %% load file and gather parameters
    load(betterLoadSave(fullpmemfilename));
    params = results;
    
    for catfile = 1:length(params)
        catfilename = extractAfter(params(catfile).filename,[bandname,'/']);
        energyfile = energyfilelist(contains({energyfilelist.name},catfilename));
        fullenergyfilename = fullfile(energyfile.folder,energyfile.name);
        load(betterLoadSave(fullenergyfilename));
        
        h = params(catfile).parameters.h;
        J = params(catfile).parameters.J;
        nodeNum = length(h);
        vectorList = mfunc_VectorList(nodeNum);
        energies = mfunc_Energy(h,J);


        
        % Gather energy features
        params(catfile).basins = features.basins;
        params(catfile).basinInds = features.basinInds;
        params(catfile).basinCount = features.basinCount;
        params(catfile).basinSizes = features.basinSizes;
        params(catfile).basin_sd = features.basin_sd;
        params(catfile).E_diff = features.E_diff;
        params(catfile).StateNumber = features.StateNumber;
        params(catfile).BasinNumber = features.BasinNumber;
        params(catfile).Dynamics = features.Dynamics;
        params(catfile).sharedFeatures = features.sharedFeatures;

        catdiscore = strrep(strrep(extractBetween(catfilename, 'cat', '_EMU'), '-', '.'),'_','');
        params(catfile).catdiscore = str2double(catdiscore{1});
    end
    
    %% save
    subfolder = fullfile(paramfolder,extractAfter(foldername,'pMEM/'));
    if ~isfolder(subfolder)
        mkdir(subfolder);
    end
    savefilename = fullfile(subfolder,strrep(filename,'result','params'));
    save(betterLoadSave(savefilename),'params');
end
        
        
    


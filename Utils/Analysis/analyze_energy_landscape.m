function analyze_energy_landscape(fname)
% ANALYZE_ENERGY_LANDSCAPE performs ELA using the fitted pMEM data. It
% computes various important features of the energy landscape. Now using
% Ezaki code.

% Inputs:
%   fname - filename with the relevant pMEM fit information (h and J)

%% Load file + Inits
load(fname, 'results');  % Load the 'results' structure from the file
[filefolder,bandname] = fileparts(fname);
bandname = extractBefore(bandname,'_result');
savedir = fullfile(strrep(filefolder,'pMEM','Energy'),bandname);
threshold = extractAfter(filefolder,'features_');


if ~isfolder(savedir)
    mkdir(savedir);
end
Nfls = length(results);

%% Loop files and do Basic ELA
parfor rec = 1:Nfls

    info = struct;
    features = struct;

    meta = get_info(results(rec).filename);

    % Info table (as struct for now)
    info.ID = meta{1};
    info.diagnosis = meta{2};
    info.catdiscore = meta{3};
    info.epoch = meta{4};
    info.task = meta{5};
    
    % Initialize features as struct
    features.info = info; %ok
    features.basinCount = NaN; %ok
    features.basinInds = NaN; %ok
    features.basins = {}; %ok
    features.basinSizes = NaN; %ok
    features.basin_sd = NaN; %ok
    features.E_diff = NaN; %ok
    features.StateNumber = NaN; %ok
    features.BasinNumber = NaN; %ok
    features.Dynamics = NaN; %ok
    features.sharedFeatures = NaN; %ok
    %% Get the original Time Series Binary Data
    binarydata = load(betterLoadSave(strrep(strrep(results(rec).filename,'Data','BinarizedData'),'clean',[threshold,'/binarized'])));
    binarydata = binarydata.binarydata;
    binarydata = bitpack_helper(binarydata,'decode');

    %% Compute Energy
    h = results(rec).parameters.h;
    J = results(rec).parameters.J;
    E = mfunc_Energy(h, J);

    %% Calculate Local Minima
    nodeNumber = length(h);
    vectorList = mfunc_VectorList(nodeNumber);
    [LocalMinIndex, BasinGraph, AdjacentList] = mfunc_LocalMin(nodeNumber, E);
    features.basinInds = LocalMinIndex;
    Nmin = length(LocalMinIndex);
    features.basinCount = Nmin;
    close all
    basins = {};
    basinSizes = [];

    basinEnergies = E(LocalMinIndex);
    GMenergy = min(basinEnergies);
    E_diff = basinEnergies-GMenergy;
    if length(E_diff) >1
        E_diff = mean(E_diff(E_diff~=0));
    else
        E_diff = NaN;
    end
    features.E_diff = E_diff;


    for m = 1:Nmin
        % Basin
        bM = vectorList(:,LocalMinIndex(m));
        bM = (bM+1)/2;
        basins{m} = sprintf('%d', bM);

        % Basin Size
        basinSizes(m) = sum(BasinGraph(:,3) == LocalMinIndex(m));

    end
    features.basins = basins';
    features.basinSizes = basinSizes';
    features.basin_sd = std(basinSizes);



    %% Calculate Disconnectivity by Dijkstra
    [~, ~] = mfunc_DisconnectivityGraph(E, LocalMinIndex, AdjacentList);
    close all %% next iter will save these

    % State Number
    StateNumber = mfunc_GetStateNumber(binarydata);
    % Basin Number
    BasinNumber = mfunc_GetBasinNumber(StateNumber, BasinGraph, LocalMinIndex);
    % Dynamics
    Dynamics = mfunc_GetDynamics(BasinNumber, Nmin);
    features.StateNumber = StateNumber;
    features.BasinNumber = BasinNumber;
    features.Dynamics = Dynamics;

    %%
    [~,fnsv] = fileparts(results(rec).filename);
    savename = fullfile(savedir,['Energy_',fnsv]);
    features2 = struct;
    features2.features = features;
    save(savename,'-v7.3','-fromstruct',features2);


end


end


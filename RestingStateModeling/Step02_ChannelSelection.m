
%% Step02 -- Calculate Channel Entropy and SNR
%
% Computes entropy and SNR for each contact in each frequency band.
% Values are aggregated across recordings by channelEntropies.
%
% Output:
%   ChannelSelection/<task>/<pt>/entropy_based.mat
%
% Variables:
%   ent_selections.<band> : Channel | Name | Entropy | SNR
%
% Requires:
%   channelEntropies.m
%   region_names from patient Regions.mat

%% Load Region Information

regionFile = betterLoadSave( ...
    ['USER_DIR/Info/', ...
    pt, '_Regions.mat']); % file contains regions for each of the bipolar contacts

load(regionFile);

%% Directories

dataDir = betterLoadSave( ...
    ['USER_DIR/Data/', ...
    task, '/', pt, '/']);

outputDir = betterLoadSave( ...
    ['USER_DIR/ChannelSelection/', ...
    task, '/', pt, '/']);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Parallel Processing

numCores = 32;

if isempty(gcp('nocreate'))
    parpool(numCores);
end

%% Calculate Entropy and SNR

bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};

ent_selections = struct();

for b = 1:numel(bands)

    band = bands{b};

    folderPath = fullfile(dataDir, ['clean_', band]);

    ent_selections.(band) = channelEntropies( ...
        folderPath, region_names);

    fprintf('Completed %s\n', band);

end

%% Save

methodstr = 'entropy';

save(fullfile(outputDir, 'entropy_based.mat'), ...
    'ent_selections', 'methodstr');

%% Generate ROI Reconstruction Channel Metadata

chanNames = { ...
    'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC', ...
    'l-dACC','l-MTG', ...
    'r-vmPFC','r-dmPFC','r-vlPFC','r-dlPFC', ...
    'r-dACC','r-MTG'};

% Preserve original band ordering
reconBands = {'alpha','beta','delta','full', ...
              'gamma','higamma','raw','theta'};

% ROI channel indices and names
recon_selections_table = table( ...
    (1:numel(chanNames))', ...
    chanNames', ...
    'VariableNames', {'Channel','Name'});

% Generate selections for each band
final_selections = struct();

for b = 1:numel(reconBands)

    final_selections(b).band = reconBands{b};
    final_selections(b).channels = recon_selections_table;

end

% Save
save(fullfile(outputDir, 'final_recon.mat'), ...
    'final_selections');

%% Close Parallel Pool

delete(gcp('nocreate'));

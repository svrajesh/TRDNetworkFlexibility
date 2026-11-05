
%% Step02 -- Weighted ROI Activity Reconstruction
%
% Reconstructs one power time series per ROI using contact-level
% entropy and SNR weights.
%
% For each band:
%   1. Load channel entropy and SNR measurements.
%   2. Load each recording's contact-level power time series.
%   3. Reconstruct ROI signals using reconstructROIweighted.
%   4. Save reconstructed signals to WeightedReconData.
%
% Output:
%   WeightedReconData/<task>/<pt>/clean_<band>/*.mat
%
% Saved variables:
%   data      : nTime x nROI reconstructed power
%   timevec   : original time vector
%   chanNames : ROI names

%% Settings

bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};

chanNames = { ...
    'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC', ...
    'l-dACC','l-MTG', ...
    'r-vmPFC','r-dmPFC','r-vlPFC','r-dlPFC', ...
    'r-dACC','r-MTG'};

%% Load Entropy and Region Information

selectionFile = betterLoadSave( ...
    ['USER_DIR/ChannelSelection/', ...
    task, '/', pt, '/entropy_based.mat']);

load(selectionFile, 'ent_selections');

regionFile = betterLoadSave( ...
    ['USER_DIR/Info/', ...
    pt, '_Regions.mat']);

load(regionFile, 'region_labels');

roilist = unique(region_labels);
roilist = roilist(roilist > 0);

%% Data Directory

datafolder = betterLoadSave( ...
    ['USER_DIR/Data/', ...
    task, '/', pt]);

%% Reconstruct Each Frequency Band

for b = 1:numel(bands)

    band = bands{b};

    band_data_folder = fullfile(datafolder, ['clean_', band]);

    band_data_files = dir(fullfile(band_data_folder, '*.mat'));
    band_data_files = band_data_files( ...
        ~startsWith({band_data_files.name}, '.'));

    % Entropy and SNR for this band
    filter_ent = ent_selections.(band);

    % Output directory
    wtRecon_data_folder = strrep( ...
        band_data_folder, '/Data/', '/WeightedReconData/');

    if ~exist(wtRecon_data_folder, 'dir')
        mkdir(wtRecon_data_folder);
    end

    %% Reconstruct Each Recording

    for f = 1:numel(band_data_files)

        filename = band_data_files(f).name;
        filePath = fullfile(band_data_folder, filename);

        load(filePath, 'data', 'timevec');

        % Entropy/SNR-weighted ROI reconstruction
        data = reconstructROIweighted( ...
            data, filter_ent, region_labels, roilist);

        % Save reconstructed ROI power time series
        save(fullfile(wtRecon_data_folder, filename), ...
            'data', 'timevec', 'chanNames', '-v7.3');

        clear data;

    end

    fprintf('Completed %s\n', band);

end

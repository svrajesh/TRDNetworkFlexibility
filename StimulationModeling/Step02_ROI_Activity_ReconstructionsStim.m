%% Step02 -- Reconstruct ROI activity for stimulation recordings
% Uses ROI weights estimated from resting-state data.
% Step01 has already segmented and labeled stimulation ON/OFF windows.

root = 'USER_DIR';

selectionFile = betterLoadSave(fullfile(root, 'ChannelSelection', 'Resting', pt, 'entropy_based.mat'));
regionFile    = betterLoadSave(fullfile(root, 'Info', [pt, '_Regions.mat']));

load(selectionFile, 'ent_selections');
load(regionFile, 'region_labels');

bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};

% Match region IDs to their original anatomical labels.
allRegionNames = {'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC','l-dACC','l-MTL', ...
                  'r-vmPFC','r-dmPFC','r-vlPFC','r-dlPFC','r-dACC','r-MTL'};
roilist = unique(region_labels(region_labels > 0));
roilist = roilist(:)';
assert(all(ismember(roilist, 1:numel(allRegionNames))), ...
    'ROI labels do not match the expected 1-12 region scheme.');
chanNames = allRegionNames(roilist);

for b = 1:numel(bands)
    band = bands{b};

    inputFolder = betterLoadSave(fullfile(root, 'Data', task, pt, ['clean_', band]));
    outputFolder = betterLoadSave(fullfile(root, 'WeightedReconData', task, pt, ['clean_', band]));
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    files = dir(fullfile(inputFolder, '*.mat'));
    files = files(~startsWith({files.name}, '.'));
    filter_ent = ent_selections.(band);

    for f = 1:numel(files)
        inputFile = fullfile(files(f).folder, files(f).name);
        S = load(inputFile, 'data', 'timevec');

        % Reconstruct all available ROIs using resting-derived weights.
        data = reconstructROIweighted(S.data, filter_ent, region_labels, roilist);
        timevec = S.timevec;

        assert(size(data, 2) == numel(chanNames), ...
            'ROI count and output channel labels do not match: %s', inputFile);

        outputFile = fullfile(outputFolder, files(f).name);
        save(outputFile, 'data', 'timevec', 'chanNames', '-v7.3');
    end
end

%% Step01 -- Preprocess Raw Neural Data
% Initialize variable pt external to this script. Ex: pt = 'DBSTRDXXX';
% Parallel Enabled

% Set up parallel processing
numCores = 32;
if isempty(gcp('nocreate'))
    parpool(numCores);
end

%% Initialize Data Directories
% Location of raw data
% Assuming raw files are .mat format
filesNS3 = 'path to NS3 data'%
filesNS5 = 'path to NS5 data'%
filesNS3 = filesNS3(~startsWith({filesNS3.name},'.'));
filesNS5 = filesNS5(~startsWith({filesNS5.name},'.'));

% Extract resting state recordings
filesRestingNS3 = filesNS3(contains({filesNS3.name},{'EyesOpen','Baseline'}));
filesRestingNS5 = filesNS5(contains({filesNS5.name},{'EyesOpen','Baseline'}));
filesCATDINS3 = filesNS3(contains({filesNS3.name},{'CATDI','CAT-DI'}));
filesCATDINS5 = filesNS5(contains({filesNS5.name},{'CATDI','CAT-DI'}));

indskeep = [];

% If both ns3 and ns5 exist, choose ns3
for i = 1:length(filesRestingNS5)
    nm = strrep(filesRestingNS5(i).name,'.ns5','.ns3');
    if any(ismember({filesRestingNS3.name},nm))
        continue
    end
    indskeep = [indskeep i];
end
filesResting = [filesRestingNS3; filesRestingNS5(indskeep)];

indskeep = [];
for i = 1:length(filesCATDINS5)
    nm = strrep(filesCATDINS5(i).name,'.ns5','.ns3');
    if any(ismember({filesCATDINS3.name},nm))
        continue
    end
    indskeep = [indskeep i];
end
filesCATDI = [filesCATDINS3; filesCATDINS5(indskeep)];

% sort files
[~,order] = sort({filesResting.name});
filesResting = filesResting(order);
filesResting = filesResting(~startsWith({filesResting.name},'.'));

% sort files
[~,order] = sort({filesCATDI.name});
filesCATDI = filesCATDI(order);
filesCATDI = filesCATDI(~startsWith({filesCATDI.name},'.'));


% Define frequency bands and initialize results structure
bands = {'raw','full','delta','theta','alpha','beta','gamma','higamma'};
files = [filesResting; filesCATDI];
files = files(~contains({files.name},'TEST'));
results = cell(1, length(files));

%% Main Loop: Loop through each file
parfor i = 1:length(files)

    % Extract filename info and generate the name for saving datafile
    filePath = fullfile(files(i).folder, files(i).name);
    info = get_info(filePath); % get_info.m returns file info as a cell array
    ID = info{1};
    diagnosis = info{2}; % TRD for now but open for other dx
    catdiscore = strrep(sprintf('%.1f', info{3}), '.', '-');
    epoch = info{4};
    taskname = info{5};
    if isempty(catdiscore)
        continue
    end
    if isempty(taskname)
        disp(info);
    end
    new_name = strcat(ID, '_', diagnosis, '_cat', catdiscore, '_', epoch,'_',taskname);



    %% Main Preprocessing Step
    dataStruct = preprocess2(filePath, 'all');  % struct with all bands


    %% Collect results for this file across all bands
    % This step is memory inefficient and may require, for large data
    % sizes, a larger RAM when parallel processing.
    fileResults = cell(1, numel(bands));
    for b = 1:numel(bands)
        %disp(strcat(num2str(i),' ',bands{b}))
        band = bands{b};
        if contains(files(i).name,{'CATDI','CAT-DI'})
            outputFolder = betterLoadSave(fullfile('USER_DIR/EnergyLandscapes/Data/CATDI', pt, ['clean_' band]));
        elseif contains(files(i).name,{'EyesOpen','Baseline'})
            outputFolder = betterLoadSave(fullfile('USER_DIR/EnergyLandscapes/Data/Resting', pt, ['clean_' band]));
        end
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end

        % Create new filename
        save_to = fullfile(outputFolder, [new_name, '.mat']);

        % Store for later saving
        fileResults{b} = {save_to, dataStruct.(band)};
    end

    results{i} = fileResults;
end

%% Save Data
% Save outside parfor
for i = 1:length(files)
    fileResults = results{i};
    if isempty(fileResults)
        continue
    end
    for b = 1:numel(bands)
        rescell = fileResults{b};

        if ~isempty(rescell)
            data = rescell{2};
            chanNames = data.Properties.VariableNames;
            data = table2array(data);
            timevec = data(:,1);
            data = data(:,2:end);
            save(rescell{1}, 'data','timevec','chanNames','-v7.3');
        end
    end
end
% Close parallel pool
delete(gcp); 




%% Step01 -- Preprocess Raw Neural Data for Stim Expt Amp Mod
% Initialize variable pt external to this script. 


%% Initialize Data Directories
% Location of raw data
% Assuming raw files are .mat format
ampmodfiles = betterLoadSave(['USER_DIR/Info/',task,'Files.xlsx']); %excel list with all the files for amp mod experiments
ptfile = readtable(ampmodfiles,'Sheet',pt);
files = cellfun(@betterLoadSave,ptfile.Filename,'UniformOutput',false);

% Define frequency bands and initialize results structure
bands = {'raw','full','delta','theta','alpha','beta','gamma','higamma'};

%% Main Loop: Loop through each file
for i = 1:length(files)
    
    % Extract filename info and generate the name for saving datafile
    filePath = files{i};%fullfile(files(i).folder, files(i).name);
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
    
    
    %% Collect results for this file across all bands
    % Since the stim recs are longer we wont par for and will save per band
    %fileResults = cell(1, numel(bands));
    for b = 1:numel(bands)
        band = bands{b};
        if contains(files{i},{'Amp'})
            outputFolder = betterLoadSave(fullfile('USER_DIR/Data/',task, pt, ['clean_' band]));
        end
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        
        % Create new filename
        save_to = fullfile(outputFolder, [new_name, '.mat']);
        [dataStruct,analysiswindows] = preprocess2(filePath, band,1);  % struct with all bands
        
        
        chanNames = dataStruct.Properties.VariableNames;
        data = table2array(dataStruct);
        timevecfull = data(:,1);
        datafull = data(:,2:end);
        
        %%
        % Parameters
        trimSamples   = 3 * 2000;       % 3 seconds trimmed off each side
        segmentLength = 60 * 2000;      % 1-minute segments (in samples)
        analysiswindowsds = ceil(analysiswindows/15);
        
        segmentwindowsds = []; % initialize output
        
        
        for w = 1:size(analysiswindowsds, 1)
            % Get trimmed window bounds
            startIdx = analysiswindowsds(w,1) + trimSamples;
            endIdx   = analysiswindowsds(w,2) - trimSamples;
            
            % Total samples in trimmed window
            L = endIdx - startIdx + 1;
            
            % How many ≥1-minute segments?
            nSeg = floor(L / segmentLength);
            
            % If window < 1 minute, keep as a single segment
            if nSeg < 1
                segmentwindowsds = [segmentwindowsds; startIdx endIdx w];
                continue
            end
            
            % Evenly spaced boundaries, INCLUSIVE indices
            rawBounds = round(linspace(startIdx, endIdx, nSeg + 1));
            
            % Convert raw boundaries to non-overlapping segments:
            %
            % Segment i:
            %   start = rawBounds(i)
            %   end   = rawBounds(i+1) - 1
            % Last segment’s end must be endIdx exactly.
            %
            segmentStarts = rawBounds(1:end-1);
            segmentEnds   = rawBounds(2:end) - 1;
            segmentEnds(end) = endIdx;   % ensure last ends exactly at endIdx
            
            % Combine into [start end parentWindowIndex]
            segmentPairs = [segmentStarts(:) segmentEnds(:) repmat(w, nSeg, 1)];
            
            % Append to output
            segmentwindowsds = [segmentwindowsds; segmentPairs];
        end
        
        
        for w = 1:size(segmentwindowsds,1)
            % Extract the segment data for the current window
            stimwin = mod(segmentwindowsds(w,3),2);
            
            switch stimwin
                case 0
                    wintype = 'ON';
                case 1
                    wintype = 'OFF';
            end
            
            data = datafull(segmentwindowsds(w, 1):segmentwindowsds(w, 2), :);
            timevec = timevecfull(segmentwindowsds(w, 1):segmentwindowsds(w, 2));
            
            save_to_final = strrep(save_to,'.mat',['_AnalysisWindow_',num2str(w),wintype,'.mat']);
            save(save_to_final, 'data','timevec','chanNames','-v7.3');
            
        end
        
        clearvars data dataStruct
        
    end
end


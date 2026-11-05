function pMEM_model_fits(base,pt,task,threshold,selections,method,ablation)

% PMEM_MODEL_FITS is the first wrapper for fitting PMEM. It passes folder
% and window information to a subfunction RUN_PMEM, and saves the results from the
% subfunction to the appropriate files

% Inputs:
%   base - root path of the project
%   pt - ID of patient
%   threshold - relevant binarization threshold %%window in MILLISECONDS
%   selections - table of the contacts selected for analysis in Step_02

% Parallel Enabled

% Set up the output directories
switch method
    % prev had other methods, but now only doing weighted reconstruction
    case 'weighted'
        data_directory = fullfile(base,'WeightedReconData',task,pt);

end
if ablation == 0
    analysis_directory = fullfile(base, 'Analysis',task,pt);
else
    analysis_directory = fullfile(base, 'Analysis','Ablation',task,pt);
end
if ablation == 0
    fold_pMEM = fullfile(analysis_directory, 'pMEM');
    fold_pMEM_feat = fullfile(fold_pMEM, ['pMEM_features_',num2str(round(threshold))]);
else
    fold_pMEM = fullfile(analysis_directory, 'pMEM');
    fold_pMEM_feat = fullfile(fold_pMEM, ['reg_',num2str(ablation),'_pMEM_features_',num2str(round(threshold))]);
end
if ~isfolder(base)
    mkdir(base);
end

if ~isfolder(fold_pMEM)
    mkdir(fold_pMEM);
end
if ~isfolder(fold_pMEM_feat)
    mkdir(fold_pMEM_feat);
end

bands = {'alpha', 'beta', 'delta', 'full', 'gamma', 'higamma', 'raw', 'theta'};


%% Loop through bands
for b = 1:length(bands)
    bname = bands{b};
    band_path = fullfile(data_directory,['clean_',bname]);
    data_files = dir(fullfile(band_path,'*.mat'));
    data_files = data_files(~startsWith({data_files.name},'.'));
    nFiles = length(data_files);

    % Process the data for each frequency band
    thisbandselections = selections(find(strcmp({selections.band},bname)));
    thisbandnames = thisbandselections.channels.Name;
    thisbandnames = cellstr(strrep(thisbandnames,"'",""));
    rfields = {'filename','parameters','metrics'};
    c = cell(length(rfields),1);
    results = cell2struct(c,rfields);

    % Loop over the data files and process each one
    parfor i = 1:nFiles
        % Load data
        data_table = load(fullfile(data_files(i).folder, data_files(i).name));  % Load the table
        chanNames = data_table.chanNames;
        if contains(chanNames,'Time')
            chanNames = chanNames(2:end);
        end
        data_table = data_table.data;
        [~,chanInds] = ismember(thisbandnames,chanNames);
        data_signals = data_table(:, chanInds);  % Extracts specified columns

        % Run pMEM analysis
        trunc = 1;
        results(i) = run_pMEM(data_signals, threshold,fullfile(data_files(i).folder, data_files(i).name),trunc,ablation);

    end

    % Save the results
    save_path = fullfile(fold_pMEM_feat, [bname, '_result.mat']);
    save(save_path,'results','-v7.3');
end
end
%% Predictive Mutual Info/Active Info Storage Runner
patients = {};
bands = {'raw'}; % run on raw data focusing on the same regions that we used before
keepreg = [1,3,4,5,7,9,10,11]; %[1,2,4,5,6,8,9,11,12,13];

for p = 1:length(patients)

    pat = patients{p};
    reginfo = load(betterLoadSave(['USER_DIR/Info/',pat,'_Regions.mat']));
    for b = 1:length(bands)

        band = bands{b};

        datafolder = betterdir(betterLoadSave(['USER_DIR/Data/Resting/',pat,'/clean_',band,'/']));


        % Preallocate results container
        results = cell(length(datafolder),1);
        regnames = reginfo.region_names;
        for rq = 1:length(regnames)
            regnames{rq} = char(strrep(regnames{rq},"'",""));
        end
        indskp = find(ismember(reginfo.region_labels,keepreg));
        regnames = regnames(indskp);


        parfor fi = 1:length(datafolder)

            fn = fullfile(datafolder(fi).folder,...
                          datafolder(fi).name);

            cat = get_cat(fn);

            S = load(fn);

            % Change "data" below to your actual variable name
            data = S.data';

            datachans = S.chanNames(2:end)';
            keepelecs = find(ismember(datachans,regnames));

            res = predictiveInformationHorizon( ...
                data(keepelecs,:),...
                2000,...
                200,...
                40,...   % history up to 200 ms
                3,...   % predict up to 15 ms ahead
                10);


            res.cat = cat;

            results{fi} = res;

        end


        % Convert cell array to structure array

        resultsAll = [results{:}];


        % Add metadata

        for fi = 1:length(resultsAll)

            resultsAll(fi).patient = pat;
            resultsAll(fi).band = band;

        end


        % Save

        svdir = betterLoadSave(fullfile('USER_DIR/Analysis/Predictability/Resting',pat));

        if ~isfolder(svdir)
            mkdir(svdir);
        end


        save(fullfile(svdir,...
            ['NewEntropy_', band, '.mat']),...
            "resultsAll",...
            '-v7.3');


        clear results resultsAll

    end

end
return


%% Extract summary metrics

n = length(resultsAll);

cat = nan(n,1);

rawPI = nan(n,1);
normPI = nan(n,1);
histTime = nan(n,1);
futureTime = nan(n,1);
areaPI = nan(n,1);

for i = 1:n

    cat(i) = resultsAll(i).cat;

    rawPI(i) = resultsAll(i).rawPI_max;

    normPI(i) = resultsAll(i).normPI_max;

    histTime(i) = resultsAll(i).historyTimescale_ms;

    futureTime(i) = resultsAll(i).futureTimescale_ms;

    areaPI(i) = resultsAll(i).rawPI_area;

end


%% Plot summary metrics

figure
tiledlayout(2,3,'TileSpacing','compact')


metrics = {
    rawPI,      'Raw predictive information (bits)';
    normPI,     'Normalized predictive information';
    areaPI,     'Predictive information area';
    histTime,   'History timescale (ms)';
    futureTime, 'Future prediction horizon (ms)'
};


for i = 1:size(metrics,1)

    nexttile

    y = metrics{i,1};

    scatter(cat,y,60,'filled')

    lsline

    xlabel('Severity')
    ylabel(metrics{i,2})

    grid on

end

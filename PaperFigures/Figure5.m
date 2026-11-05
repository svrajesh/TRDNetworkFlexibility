%% Figure 5
%% Gather Info
allpts = {};
alllabels = {};
alltasks = {'Resting'};
allbands = {'delta','theta','alpha','beta','gamma'};
allregions = {'l-vmPFC','l-vlPFC','l-dlPFC','l-dACC', ...
              'r-vmPFC','r-vlPFC','r-dlPFC','r-dACC'};
thresholds = [50];
USER_DIR = '';

outDir = fullfile(rootDir,'PaperFigures');

ptgroups        = [];
bandgroups      = [];
reggroups       = [];
recordinggroups = [];

aperiodicsAll = [];
tempsAll      = [];
R2pmemAll     = [];
scoresAll     = [];
transPall     = [];
dtAll         = [];
entAll        = [];

for pt = 1:length(allpts)

    for tk = 1:length(alltasks)

        for b = 1:length(allbands)

            for thresh = 1:length(thresholds)

                paramfile = betterLoadSave([USER_DIR,'/Analysis/',alltasks{tk},'/',allpts{pt},'/Params/pMEM_features_',num2str(thresholds(thresh)),'/',allbands{b},'_params.mat']);
                load(paramfile);

                %% ---------------------------------------------
                % Recording index
                %
                % pf indexes the recording within this
                % patient/band.
                %% ---------------------------------------------
                recordinggrouping = (1:length(params))';

                %% ---------------------------------------------
                % Loop over regions
                %% ---------------------------------------------
                for r = 1:length(allregions)

                    ptgrouping  = pt * ones(length(params),1);
                    bandgrouping = b * ones(length(params),1);
                    reggrouping = r * ones(length(params),1);

                    %% -----------------------------------------
                    % Extract parameters
                    %% -----------------------------------------
                    thisApers = [params.parameters];

                    thisTemps = [thisApers.effectiveTemperature];

                    thisApers = [thisApers.aperiodicSlopes];

                    thisApers = thisApers(r,:);

                    %% -----------------------------------------
                    % Extract PMEM metrics
                    %% -----------------------------------------
                    mets = [params.metrics];
                    mets = [mets.sharedMetrics];

                    thisAcc = [mets.d];
                    thisR2pmem = [mets.R2pmem];

                    %% -----------------------------------------
                    % Transition dynamics
                    %% -----------------------------------------
                    Rall = [];
                    meanDwellAll = [];
                    HrateAll = [];
                    NeffAll = [];

                    for pf = 1:length(params)

                        Dtrans = ...
                            params(pf).sharedFeatures.Dynamics.Dtrans;

                        [ ...
                            Rall(pf), ...
                            meanDwellAll(pf), ...
                            HrateAll(pf), ...
                            NeffAll(pf)] = ...
                            summarizeTransitions(Dtrans);

                    end

                    thisR  = Rall;
                    thisDT = meanDwellAll;
                    thisEnt = HrateAll;

                    %% -----------------------------------------
                    % Symptom scores
                    %% -----------------------------------------
                    thisScores = [params.catdiscore];

                    %% -----------------------------------------
                    % Append grouping variables
                    %% -----------------------------------------
                    ptgroups = [ ...
                        ptgroups; ...
                        ptgrouping(:)];

                    bandgroups = [ ...
                        bandgroups; ...
                        bandgrouping(:)];

                    reggroups = [ ...
                        reggroups; ...
                        reggrouping(:)];

                    recordinggroups = [ ...
                        recordinggroups; ...
                        recordinggrouping(:)];

                    %% -----------------------------------------
                    % Append features
                    %% -----------------------------------------
                    aperiodicsAll = [ ...
                        aperiodicsAll; ...
                        thisApers(:)];

                    tempsAll = [ ...
                        tempsAll; ...
                        thisTemps(:)];

                    R2pmemAll = [ ...
                        R2pmemAll; ...
                        thisR2pmem(:)];

                    scoresAll = [ ...
                        scoresAll; ...
                        thisScores(:)];

                    dtAll = [ ...
                        dtAll; ...
                        thisDT(:)];

                    transPall = [ ...
                        transPall; ...
                        thisR(:)];

                    entAll = [ ...
                        entAll; ...
                        thisEnt(:)];

                end
            end
        end
    end
end
statsDir = [USER_DIR,'StatTables'];

exportFigure5Statistics(aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,reggroups,recordinggroups,R2pmemAll, ...
    allpts,allregions,alllabels,fullfile(statsDir,'Figure5_Statistics.xlsx'));

plotMetricFig5FullNew( ...
    aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,reggroups,recordinggroups,R2pmemAll, ...
    allpts,allregions,alllabels,3.5)

exportgraphics( ...
    gcf, ...
    betterLoadSave(fullfile(outDir,'Fig5_ScatterPlots.svg')), ...
    'ContentType','vector');
plotMetricFig5FullNew( ...
    aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,reggroups,recordinggroups,R2pmemAll, ...
    allpts,allregions,alllabels,7.5)

exportgraphics( ...
    gcf, ...
    betterLoadSave(fullfile(outDir,'Fig5_HeatmapsForests.svg')), ...
    'ContentType','vector');
close all

allpts = {};
alllabels = {};
alltasks = {'Resting'};
allbands = {'delta','theta','alpha','beta','gamma'};
thresholds = 50;
USER_DIR = '';
outDir = fullfile(USER_DIR,'PaperFigures');
statsDir = fullfile(USER_DIR,'StatTables');
if ~isfolder(statsDir)
    mkdir(statsDir);
end

ptgroups = [];
recordinggroups = [];
bandgroups = [];

accuracyAll = [];
R2pmemAll = [];
tempsAll = [];
scoresAll = [];
transPall = [];
dtAll = [];
entAll = [];

recordingCounter = 0;

for pt = 1:length(allpts)
    for tk = 1:length(alltasks)
        for b = 1:length(allbands)
            for thresh = 1:length(thresholds)

                paramfile = betterLoadSave([USER_DIR,'/Analysis/',alltasks{tk},'/',allpts{pt},'/Params/pMEM_features_',num2str(thresholds(thresh)),'/',allbands{b},'_params.mat']);
                load(paramfile);

                nRecordings = length(params);

                % Globally unique recording IDs.
                %recordingIDs = recordingCounter + (1:nRecordings)';
                recordingIDs = (1:nRecordings)';
                recordingCounter = recordingCounter + nRecordings;

                ptgrouping = pt * ones(nRecordings,1);
                bandgrouping = b * ones(nRecordings,1);

                % ---------------------------------------------------------
                % Extract metrics
                % ---------------------------------------------------------
                mets = [params.metrics];
                mets = [mets.sharedMetrics];

                thisAcc = [mets.d];
                thisR2pmem = [mets.R2pmem];

                thisTemps = [params.parameters];
                thisTemps = [thisTemps.effectiveTemperature];

                Rall = nan(nRecordings,1);
                meanDwellAll = nan(nRecordings,1);
                HrateAll = nan(nRecordings,1);
                NeffAll = nan(nRecordings,1);

                for pf = 1:nRecordings
                    Dtrans = params(pf).sharedFeatures.Dynamics.Dtrans;

                    [Rall(pf), ...
                     meanDwellAll(pf), ...
                     HrateAll(pf), ...
                     NeffAll(pf)] = summarizeTransitions(Dtrans);
                end

                thisR = Rall;
                thisDT = meanDwellAll;
                thisEnt = HrateAll;
                thisScores = [params.catdiscore];

                % ---------------------------------------------------------
                % Append
                % ---------------------------------------------------------
                ptgroups = [ptgroups; ptgrouping];
                recordinggroups = [recordinggroups; recordingIDs];
                bandgroups = [bandgroups; bandgrouping];

                accuracyAll = [accuracyAll; thisAcc(:)];
                R2pmemAll = [R2pmemAll; thisR2pmem(:)];

                scoresAll = [scoresAll; thisScores(:)];
                tempsAll = [tempsAll; thisTemps(:)];
                dtAll = [dtAll; thisDT(:)];
                transPall = [transPall; thisR(:)];
                entAll = [entAll; thisEnt(:)];

            end
        end
    end
end
exportFigure3Statistics(scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels, ...
    fullfile(statsDir,'Figure3_Statistics.xlsx'));

q1 = plotMetricFig3FullNew( ...
    scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels,3.5);
exportgraphics( ...
    gcf, ...
    betterLoadSave(fullfile(outDir,'Fig3_ScatterPlots.svg')), ...
    'ContentType','vector');
q2 = plotMetricFig3FullNew( ...
    scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels,7.5);
exportgraphics( ...
    gcf, ...
    betterLoadSave(fullfile(outDir,'Fig3_HeatmapsForests.svg')), ...
    'ContentType','vector');

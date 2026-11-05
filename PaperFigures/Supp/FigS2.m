%% FIGURE S2 | CAT-DI severity versus pMEM landscape dynamics
% CATDI remains the input dataset. Matched recording IDs and a
% minimum of three R2-eligible recordings per patient/band are enforced.

allpts = {};
alllabels = {};
alltasks = {'CATDI'};
allbands = {'delta','theta','alpha','beta','gamma'};
thresholds = 50;
outDir = [USER_DIR,'/CleanFiguresNew'];

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

for pt = 1:length(allpts)
    for tk = 1:length(alltasks)
        for b = 1:length(allbands)
            for thresh = 1:length(thresholds)

                paramfile = betterLoadSave([USER_DIR,'/Analysis/',alltasks{tk},'/',allpts{pt},'/Params/pMEM_features_',num2str(thresholds(thresh)),'/',allbands{b},'_params.mat']);
                load(paramfile);

                nRecordings = length(params);

                % Match Fig. 3: the same within-patient recording index is
                % reused across frequency bands. Assumes matched ordering.
                recordingIDs = (1:nRecordings)';

                ptgrouping = pt * ones(nRecordings,1);
                bandgrouping = b * ones(nRecordings,1);

                % Extract metrics
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

                % Append
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

% Require at least three R2-eligible recordings per patient and band.
% The detailed helper already enforces this; here it also applies to the
% pooled beta-band panel, which otherwise accepts two observations.
minRecordings = 3;
R2threshold = 0.6;
keep = true(size(scoresAll));
for pt = 1:numel(allpts)
    for b = 1:numel(allbands)
        inGroup = ptgroups == pt & bandgroups == b;
        nEligible = nnz(inGroup & R2pmemAll >= R2threshold & ...
            isfinite(scoresAll));
        if nEligible < minRecordings
            keep(inGroup) = false;
        end
    end
end

ptgroups = ptgroups(keep);
recordinggroups = recordinggroups(keep);
bandgroups = bandgroups(keep);
R2pmemAll = R2pmemAll(keep);
scoresAll = scoresAll(keep);
tempsAll = tempsAll(keep);
transPall = transPall(keep);
dtAll = dtAll(keep);
entAll = entAll(keep);


statsDir = [USER_DIR,'/StatTables'];

exportFigure3Statistics(scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels, ...
    fullfile(statsDir,'FigureS2_Statistics.xlsx'));



% Both panels use the same Fig. 3 helper, with a different figure height.
q1 = plotMetricFig3FullNew( ...
    scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels,3.5);
exportgraphics( ...
    gcf, ...
    betterLoadSave(fullfile(outDir,'FigS2_ScatterPlots.svg')), ...
    'ContentType','vector');
q2 = plotMetricFig3FullNew( ...
    scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels,7.5);
exportgraphics( ...
    gcf, ...
    betterLoadSave(fullfile(outDir,'FigS2_HeatmapsForests.svg')), ...
    'ContentType','vector');

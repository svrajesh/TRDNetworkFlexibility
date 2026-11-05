%% FIGURE S1 | CAT-DI distribution, spectral correlations, and model diagnostics
% Original data selection, calculations, and inference retained.

clear

% SETUP
allpts = {};
alltasks = {'CATDI'};
alllabels = {};

chanNames = {'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC','l-dACC','l-MTG', ...
    'r-vmPFC','r-dmPFC','r-vlPFC','r-dlPFC','r-dACC','r-MTG'};

keepreg = [1,3,4,5,7,9,10,11];

bandsPower = {'delta','theta','alpha','beta','gamma','aperiodic'};
bandLabels = {'δ','θ','α','β','γ','1/f'};

allbands = {'delta','theta','alpha','beta','gamma'};
thresholds = 50;

nPts = length(allpts);
nBandsPower = length(bandsPower);
nChan = length(keepreg);

% FIGURE
fig = figure('Units','pixels', ...
    'Position',[100 50 1250 900], ...
    'Renderer','painters', ...
    'Color','w');

tl = tiledlayout(fig,4,6, ...
    'Padding','compact', ...
    'TileSpacing','tight');

% TOP LEFT: WITHIN-SUBJECT QQ PLOTS
% Rows 1:2, columns 1:3
axQQ = nexttile(tl,1,[2 3]);
hold(axQQ,'on');

cmapPts = lines(nPts);

for pt = 1:nPts

    filePath = fullfile( ...
        USER_DIR,'/Analysis/CATDI', ...
        allpts{pt}, ...
        '/Params/pMEM_features_50/alpha_params.mat');

    load(betterLoadSave(filePath),'params');

    thisScores = [params.catdiscore];
    thisScores = thisScores(:);
    thisScores = thisScores(~isnan(thisScores));

    % Standardize within subject
    thisScoresZ = zscore(thisScores);

    % QQ quantiles
    n = length(thisScoresZ);
    p = ((1:n)' - 0.5) ./ n;

    theoreticalQ = norminv(p);
    sampleQ = sort(thisScoresZ);

    scatter(axQQ, ...
        theoreticalQ, ...
        sampleQ, ...
        'filled',...
        'MarkerFaceColor',cmapPts(pt,:), ...
        'MarkerFaceAlpha',0.5,...
        'DisplayName',alllabels{pt});
end

% Reference normal line
xl = xlim(axQQ);
yl = ylim(axQQ);
lims = [min([xl yl]), max([xl yl])];

plot(axQQ,lims,lims,'k--', ...
    'LineWidth',1.5, ...
    'HandleVisibility','off');

xlim(axQQ,lims);
ylim(axQQ,lims);

axis(axQQ,'square');
box(axQQ,'off');

xlabel(axQQ,'Standard Normal Quantiles','FontSize',14);
ylabel(axQQ,'Within-Subject Standardized Scores','FontSize',14);
title(axQQ,'QQ Analysis of CAT-DI Normality','FontSize',16);

legend(axQQ, ...
    'Location','southeast', ...
    'Box','off', ...
    'FontSize',10);

set(axQQ,'FontSize',12);

% GATHER POWER CORRELATIONS
R_allPts = cell(nPts,1);
P_allPts = cell(nPts,1);

for pt = 1:nPts

    for tk = 1:length(alltasks)

        load(betterLoadSave([ ...
            USER_DIR,'/Analysis/', ...
            alltasks{tk},'/',allpts{pt}, ...
            '/PowerAnalysis/powerAnalysisWeightedRecon.mat']));

        load(betterLoadSave([ ...
            USER_DIR,'/Analysis/', ...
            alltasks{tk},'/',allpts{pt}, ...
            '/Params/pMEM_features_50/raw_params.mat']));

        R_values = zeros(nBandsPower,nChan);
        P_values = zeros(nBandsPower,nChan);

        for bandIndex = 1:nBandsPower

            band = bandsPower{bandIndex};
            scores = scoresall(bandIndex,:);

            for chanIndex = 1:nChan

                if strcmp(band,'aperiodic')

                    apers = [params.parameters];
                    apers = [apers.aperiodicSlopes];

                    bandData = apers(chanIndex,:);

                else

                    bandpowerbc = squeeze( ...
                        bandpower(bandIndex,:,keepreg(chanIndex)));

                    bandData = 20*log10(bandpowerbc);

                end

                [R,P] = corrcoef(bandData,scores);

                R_values(bandIndex,chanIndex) = R(1,2);
                P_values(bandIndex,chanIndex) = P(1,2);

            end
        end

        R_allPts{pt} = R_values;
        P_allPts{pt} = P_values;

    end
end

% TOP RIGHT: POWER CORRELATION HEATMAPS
% Rows 1:2, columns 4:6
for pt = 1:nPts

    % 2 x 3 arrangement
    row = ceil(pt/3);
    col = mod(pt-1,3) + 4;

    tileNum = (row-1)*6 + col;

    ax = nexttile(tl,tileNum);
    hold(ax,'on');

    R_values = R_allPts{pt};
    P_values = P_allPts{pt};

    nRows = size(R_values,1);
    nCols = size(R_values,2);

    % pcolor logic matching the working code
    R_plot = nan(nRows+1,nCols+1);
    R_plot(1:nRows,1:nCols) = R_values;

    xEdges = 0.5:1:(nCols+0.5);
    yEdges = 0.5:1:(nRows+0.5);

    pcolor(ax,xEdges,yEdges,R_plot);
    shading(ax,'flat');

    set(ax,'YDir','reverse');

    colormap(ax,redblue);

    ax.CLim = [-1 1];

    xlim(ax,[0.5 nCols+0.5]);
    ylim(ax,[0.5 nRows+0.5]);

    title(ax,alllabels{pt},'FontSize',11);

    if row == 2
        xticklabels(ax,chanNames(keepreg));
        xtickangle(ax,45);
        %% X axis
        set(ax, ...
            'XTick',1:nCols, ...
            'FontSize',10);

    else
        xticklabels(ax,[]);
    end

    %% Y axis
    set(ax,'YTick',1:nRows);

    if col == 4
        ylabel(ax,'Bands');
        yticklabels(ax,bandLabels);
    else
        yticklabels(ax,[]);
    end

    box(ax,'off');

    %% FDR correction
    fdr = fdr_BH(reshape(P_values,1,[]),0.05);
    pfdr = reshape(fdr,size(P_values));

    for b = 1:size(pfdr,1)

        for r = 1:size(pfdr,2)

            if pfdr(b,r) < 0.05

                text(ax,r,b,'*', ...
                    'Color','k', ...
                    'FontSize',11, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle');

            end
        end
    end

end

% % Shared colorbar
% axCB = nexttile(tl,6);
% cb = colorbar(axCB);
% cb.Label.String = 'Pearson r';
% cb.Label.FontSize = 11;

% GATHER MODEL EVALUATION DATA
ptgroups = [];
bandgroups = [];
accuracyAll = [];
R2pmemAll = [];

tempsAll = [];
scoresAll = [];
transPall = [];
dtAll = [];

for pt = 1:nPts

    for tk = 1:length(alltasks)

        for b = 1:length(allbands)

            for thresh = 1:length(thresholds)

                paramfile = betterLoadSave([ ...
                    USER_DIR,'/Analysis/', ...
                    alltasks{tk},'/',allpts{pt}, ...
                    '/Params/pMEM_features_',num2str(round(thresholds(thresh))),'/', ...
                    allbands{b},'_params.mat']);

                load(paramfile);

                ptgrouping = pt*ones(length(params),1);
                bandgrouping = b*ones(length(params),1);

                mets = [params.metrics];
                mets = [mets.sharedMetrics];

                thisAcc = [mets.d];
                thisR2pmem = [mets.R2pmem];

                thisTemps = [params.parameters];
                thisTemps = [thisTemps.effectiveTemperature];

                Rall = [];
                meanDwellAll = [];
                HrateAll = [];
                NeffAll = [];

                for pf = 1:length(params)

                    Dtrans = params(pf).sharedFeatures.Dynamics.Dtrans;

                    [Rall(pf), ...
                        meanDwellAll(pf), ...
                        HrateAll(pf), ...
                        NeffAll(pf)] = summarizeTransitions(Dtrans);

                end

                thisR = Rall;
                thisDT = meanDwellAll;
                thisScores = [params.catdiscore];

                ptgroups = [ptgroups; ptgrouping(:)];
                bandgroups = [bandgroups; bandgrouping(:)];

                accuracyAll = [accuracyAll; thisAcc(:)];
                R2pmemAll = [R2pmemAll; thisR2pmem(:)];

                scoresAll = [scoresAll; thisScores(:)];
                tempsAll = [tempsAll; thisTemps(:)];
                dtAll = [dtAll; thisDT(:)];
                transPall = [transPall; thisR(:)];

            end
        end
    end
end

rgbMatrix = cmapPts(ptgroups,:);

% BOTTOM LEFT: R2
% ONLY COLUMN 1
axR2 = nexttile(tl,13,[2,1]);
hold(axR2,'on');

violinplot(axR2, ...
    bandgroups, ...
    R2pmemAll, ...
    "FaceColor","none", ...
    "EdgeColor","k", ...
    "LineWidth",1.5);

swarmchart(axR2, ...
    bandgroups, ...
    R2pmemAll, ...
    7, ...
    rgbMatrix, ...
    'filled', ...
    'MarkerFaceAlpha',0.5);

yline(axR2,0.6,'--', ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5);

title(axR2,'Model Fit','FontSize',14);
ylabel(axR2,'R^2','FontSize',13);

xticks(axR2,1:5);
xticklabels(axR2,{'δ','θ','α','β','γ'});

ylim(axR2,[0 1]);
box(axR2,'off');
set(axR2,'FontSize',11);

% BOTTOM: ACCURACY
% ONLY COLUMN 2
axAcc = nexttile(tl,14,[2,1]);
hold(axAcc,'on');

violinplot(axAcc, ...
    bandgroups, ...
    accuracyAll, ...
    "FaceColor","none", ...
    "EdgeColor","k", ...
    "LineWidth",1.5);

swarmchart(axAcc, ...
    bandgroups, ...
    accuracyAll, ...
    12, ...
    rgbMatrix, ...
    'filled', ...
    'MarkerFaceAlpha',0.7);

title(axAcc,'Model Accuracy','FontSize',14);
ylabel(axAcc,'r_D','FontSize',13);

xticks(axAcc,1:5);
xticklabels(axAcc,{'δ','θ','α','β','γ'});

ylim(axAcc,[0 1]);
box(axAcc,'off');
set(axAcc,'FontSize',11);

% THRESHOLD ANALYSIS
% Uses the 4 x 6 tiledlayout created above.

allpts = {'DBSTRD002','DBSTRD006','DBSTRD008','DBSTRD010','DBSTRD011','DBSTRD014'};
alllabels = {'Patient 1','Patient 2','Patient 3','Patient 4','Patient 5','Patient 6'};
alltasks = {'Resting'};

allbands = {'delta','theta','alpha','beta','gamma'};
allthresholds = [25,30,35,37.5,40,45,50,55,60,62.5,65,70,75];

% Gather data

ptgroups = [];
bandgroups = [];
threshgroups = [];
accuracyAll = [];
R2pmemAll = [];

tempsAll = [];
scoresAll = [];
transPall = [];
dtAll = [];

for pt = 1:length(allpts)

    for tk = 1:length(alltasks)

        for b = 1:length(allbands)

            for thresh = 1:length(allthresholds)

                paramfile = betterLoadSave([ ...
                    USER_DIR,'/Analysis/', ...
                    alltasks{tk},'/',allpts{pt}, ...
                    '/Params/pMEM_features_', ...
                    num2str(round(allthresholds(thresh))), ...
                    '/',allbands{b},'_params.mat']);

                load(paramfile);

                % Grouping variables
                ptgrouping = pt * ones(length(params),1);
                bandgrouping = b * ones(length(params),1);
                threshgrouping = thresh * ones(length(params),1);

                % Model metrics
                mets = [params.metrics];
                mets = [mets.sharedMetrics];

                thisAcc = [mets.d];
                thisR2pmem = [mets.R2pmem];

                % Effective temperature
                thisTemps = [params.parameters];
                thisTemps = [thisTemps.effectiveTemperature];

                % Transition metrics
                Rall = nan(length(params),1);
                meanDwellAll = nan(length(params),1);
                HrateAll = nan(length(params),1);
                NeffAll = nan(length(params),1);

                for pf = 1:length(params)

                    Dtrans = params(pf).sharedFeatures.Dynamics.Dtrans;

                    [Rall(pf), ...
                        meanDwellAll(pf), ...
                        HrateAll(pf), ...
                        NeffAll(pf)] = summarizeTransitions(Dtrans);

                end

                thisR = Rall;
                thisDT = meanDwellAll;

                thisScores = [params.catdiscore];

                % Append
                ptgroups = [ptgroups; ptgrouping(:)];
                bandgroups = [bandgroups; bandgrouping(:)];
                threshgroups = [threshgroups; threshgrouping(:)];

                accuracyAll = [accuracyAll; thisAcc(:)];
                R2pmemAll = [R2pmemAll; thisR2pmem(:)];

                tempsAll = [tempsAll; thisTemps(:)];
                scoresAll = [scoresAll; thisScores(:)];
                dtAll = [dtAll; thisDT(:)];
                transPall = [transPall; thisR(:)];

            end
        end
    end
end

% Plot
% Last 4 columns, spanning both bottom rows

axThresh = nexttile(tl,15,[2 4]);
hold(axThresh,'on');

cmap = lines(max(ptgroups));
rgbMatrix = cmap(ptgroups,:);

rawthresholds = allthresholds;

% Average accuracy across bands assuming identical recording/threshold ordering.

threshgroupsBA = [];
accuracyAllBA = [];

for b = 1:length(allbands)

    bandinds = find(bandgroups == b);

    if b == 1

        accuracyAllBA = accuracyAll(bandinds);
        threshgroupsBA = threshgroups(bandinds);
        rgbMatrixBA = rgbMatrix(bandinds,:);

    else

        accuracyAllBA = accuracyAllBA + accuracyAll(bandinds);

    end

end

accuracyAllBA = accuracyAllBA / length(allbands);

%% Swarm + box plot

swarmchart(axThresh, ...
    threshgroupsBA, ...
    accuracyAllBA, ...
    15, ...
    rgbMatrixBA, ...
    'filled', ...
    'MarkerFaceAlpha',0.5);

boxchart(axThresh, ...
    threshgroupsBA, ...
    accuracyAllBA, ...
    'BoxFaceColor','none', ...
    'BoxEdgeColor','k', ...
    'LineWidth',2, ...
    'MarkerStyle','.');

%% Formatting

xlabel(axThresh,'Thresholds');
ylabel(axThresh,'r_D');

box(axThresh,'off');

xticks(axThresh,1:length(allthresholds));
xticklabels(axThresh,rawthresholds);

ylim(axThresh,[0 1]);

set(axThresh,'FontSize',12);

title(axThresh,'Resting State r_D Value Across Thresholds');


drawnow;
pause(0.2);
fig = gcf;
exportgraphics(fig, ...
    betterLoadSave(USER_DIR,'/PaperFigures/FigS1_Combined.svg'), ...
    'Resolution',600);

clear

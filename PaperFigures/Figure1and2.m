%% ========================================================================
% FIGURES 1 AND 2
%
% FIGURE 1
%   A. Depression severity (CAT-DI) over inpatient stay
%   B. Correlation of spectral features with CAT-DI severity
%
% FIGURE 2
%   A. Brain schematic (spans two subplot positions)
%   B. Example severity/model relationship
%   C. Energy landscape
%   D. Band-filtered signals
%   E. Binarized signals
%   F. Model R2
%   G. Model accuracy (r_D)
%
% Exports individual panels to CleanFiguresNew.
% Figure 1 exports: Fig1_*
% Figure 2 exports: Fig2_*
% ========================================================================

clear;
clc;

%% Settings

USER_DIR = '';

outDir = fullfile(rootDir,'PaperFigures');

allpts = {};

alllabels = {};

nPts = numel(allpts);
ptColors = lines(nPts);


%% ========================================================================
% FIGURE 1A: DEPRESSION SEVERITY OVER TIME
% ========================================================================

daylabels = arrayfun(@(d) sprintf('Day %d',d), ...
    1:9,'UniformOutput',false);

fig = figure('Units','pixels', ...
    'Position',[0 0 800 500], ...
    'Renderer','painters', ...
    'Color','w');

t = tiledlayout(3,2);

for pt = 1:nPts

    patientID = allpts{pt};

    % Load CAT-DI scores
    catdiFile = betterLoadSave(fullfile( ...
        rootDir,'CATDI',[patientID '_CATDI_Scores.xlsx']));

    catdi = readtable(catdiFile);

    % Keep resting-state observations
    keepinds = find(~cellfun(@isempty,catdi.Resting));

    score = catdi.Score(keepinds);
    timestamp = catdi.Time(keepinds);

    % Plot scores
    ax = nexttile(t);

    scatter(ax,timestamp,score,20,ptColors(pt,:),'filled');

    hold(ax,'on');
    xline(ax,1:8,'k:','LineWidth',1);
    hold(ax,'off');

    xlim(ax,[0 9]);
    ylim(ax,[20 80]);

    title(ax,alllabels{pt});
    box(ax,'off');

    % Only show time labels on bottom row
    if pt <= 4

        xticks(ax,[]);

    else

        xticks(ax,0.5:1:8.5);
        xticklabels(ax,daylabels);
        xlabel(ax,'Time');

        set(ax,'TickLength',[0 0]);

    end

end

ylabel(t,'Depression Severity');

sgtitle(t,'Depression Severity Scores Over Inpatient Stay', ...
    'FontSize',16);

drawnow;

exportgraphics(fig, ...
    betterLoadSave(fullfile(outDir,'Fig1_CATDI_All.svg')), ...
    'ContentType','vector');

close(fig);


%% ========================================================================
% FIGURE 1B: SPECTRAL FEATURES VS CAT-DI
%
% Six panels, one per patient.
% Correlates five frequency bands and aperiodic slope with CAT-DI.
% FDR correction is performed within each patient across 6 x 8 tests.
% ========================================================================

chanNames = {'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC', ...
             'l-dACC','l-MTG','r-vmPFC','r-dmPFC', ...
             'r-vlPFC','r-dlPFC','r-dACC','r-MTG'};

keepreg = [1 3 4 5 7 9 10 11];

bands = {'delta','theta','alpha','beta','gamma','aperiodic'};

bandLabels = {'δ','θ','α','β','γ','1/f'};

nBands = numel(bands);
nRegions = numel(keepreg);

fig = figure('Units','pixels', ...
    'Position',[0 0 800 500], ...
    'Renderer','painters', ...
    'Color','w');

t = tiledlayout(2,3, ...
    'Padding','loose', ...
    'TileSpacing','compact');

for pt = 1:nPts

    patientID = allpts{pt};

    %% Load spectral power and model parameters

    analysisDir = fullfile(rootDir,'Analysis','Resting',patientID);

    powerFile = betterLoadSave(fullfile( ...
        analysisDir,'PowerAnalysis','powerAnalysisWeightedRecon.mat'));

    paramsFile = betterLoadSave(fullfile( ...
        analysisDir,'Params','pMEM_features_50','raw_params.mat'));

    powerData = load(powerFile);
    paramData = load(paramsFile);

    bandpower = powerData.bandpower;
    scoresall = powerData.scoresall;
    params = paramData.params;

    %% Calculate correlations

    R_values = zeros(nBands,nRegions);
    P_values = zeros(nBands,nRegions);

    % Aperiodic slopes from individual recording parameters
    apers = [params.parameters];
    apers = [apers.aperiodicSlopes];

    for b = 1:nBands

        scores = scoresall(b,:);

        for r = 1:nRegions

            if strcmp(bands{b},'aperiodic')

                % Aperiodic slope data are already indexed by
                % the eight retained regions.
                bandData = apers(r,:);

            else

                % Log-transformed spectral amplitude
                bandpowerbc = squeeze(bandpower(b,:,keepreg(r)));
                bandData = 20*log10(bandpowerbc);

            end

            [R,P] = corrcoef(bandData,scores);

            R_values(b,r) = R(1,2);
            P_values(b,r) = P(1,2);

        end
    end

    %% FDR correction across all band x region comparisons

    fdr = fdr_BH(P_values(:)',0.05);
    pfdr = reshape(fdr,size(P_values));

    %% Plot heatmap

    ax = nexttile(t);

    % Pad pcolor to preserve the final row and column
    R_plot = nan(nBands+1,nRegions+1);
    R_plot(1:nBands,1:nRegions) = R_values;

    xEdges = 0.5:1:(nRegions+0.5);
    yEdges = 0.5:1:(nBands+0.5);

    pcolor(ax,xEdges,yEdges,R_plot);
    shading(ax,'flat');

    set(ax,'YDir','reverse');

    xlim(ax,[0.5 nRegions+0.5]);
    ylim(ax,[0.5 nBands+0.5]);

    clim(ax,[-1 1]);
    colormap(ax,redblue);

    title(ax,alllabels{pt});

    set(ax, ...
        'XTick',1:nRegions, ...
        'YTick',1:nBands, ...
        'FontSize',12);

    % Channel labels on bottom row
    if pt > 3

        xlabel(ax,'Channels');
        xticklabels(ax,chanNames(keepreg));

    else

        xticklabels(ax,[]);

    end

    % Band labels on left column
    if mod(pt,3) == 1

        ylabel(ax,'Bands');
        yticklabels(ax,bandLabels);

    else

        yticklabels(ax,[]);

    end

    %% Mark FDR-significant correlations

    hold(ax,'on');

    for b = 1:nBands
        for r = 1:nRegions

            if pfdr(b,r) < 0.05

                text(ax,r,b,'*', ...
                    'Color','k', ...
                    'FontSize',12, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle');

            end

        end
    end

    hold(ax,'off');

    %% Shared colorbar

    if pt == nPts

        cb = colorbar(ax);
        cb.Layout.Tile = 'east';

    end

end

sgtitle(t,'Correlation of Regional Spectral Power with CAT-DI Score', ...
    'FontSize',16);

drawnow;

exportgraphics(fig, ...
    betterLoadSave(fullfile(outDir,'Fig1_PowerCorrPlots.svg')), ...
    'ContentType','vector');

close(fig);


%% ========================================================================
% FIGURE 2: pMEM AND ENERGY LANDSCAPE SCHEMATIC
%
% The brain spans subplot positions 1-2.
% Seven distinct axes occupy the eight-position layout.
%
% Retains:
%   - Brain visualization
%   - Schematic parameter correlation
%   - Energy landscape
%   - Filtered and binarized neural signals
%   - R2 and accuracy distributions
% ========================================================================

%% Load brain and example neural signals

load(fullfile(rootDir,'Codebase','Utils','Figures','EC.mat'));
load(fullfile(rootDir,'Codebase','Utils','Figures','SURF.mat'));

signalFile = '';%example filename here

signalData = load(signalFile);

% Example recording and selected prefrontal channels
window = 67000:74999;
dat = signalData.data(window,[1 3 4 5 7 9 10 11]);

% Region colors (preserved from original figure)
cols = [
    0.5, 0.5, 0.5;
    1, 0, 0;
    .1490, .5490, .8660;
    0, 1, 0;
    1, 0, 1;
    1, 0, 0;
    .1490, .5490, .8660;
    0, 1, 0;
    1, 0, 1
];

% Normalize example signals once for both plots
datNorm = zeros(size(dat));

for i = 1:size(dat,2)

    chan = dat(:,i);
    datNorm(:,i) = (chan-min(chan))/range(chan);

end


%% ========================================================================
% Gather pMEM goodness-of-fit and accuracy
%
% Only these two metrics are required for Figure 2.
% ========================================================================

allbands = {'delta','theta','alpha','beta','gamma'};

ptgroups = [];
bandgroups = [];

R2pmemAll = [];
accuracyAll = [];

for pt = 1:nPts

    patientID = allpts{pt};

    for b = 1:numel(allbands)

        % Load pMEM parameters for each patient and frequency band
        paramfile = betterLoadSave(fullfile( ...
            rootDir,'Analysis','Resting',patientID, ...
            'Params','pMEM_features_50', ...
            [allbands{b} '_params.mat']));

        S = load(paramfile);
        params = S.params;

        nRecordings = numel(params);

        % Model fit metrics
        mets = [params.metrics];
        mets = [mets.sharedMetrics];

        thisAcc = [mets.d];
        thisR2 = [mets.R2pmem];

        % Group labels
        ptgroups = [ptgroups; repmat(pt,nRecordings,1)];
        bandgroups = [bandgroups; repmat(b,nRecordings,1)];

        % Accumulate metrics
        accuracyAll = [accuracyAll; thisAcc(:)];
        R2pmemAll = [R2pmemAll; thisR2(:)];

    end
end

% Patient colors for model performance plots
cmapPts = lines(nPts);
rgbMatrix = cmapPts(ptgroups,:);


%% ========================================================================
% Create Figure 2
% ========================================================================

fig = figure( ...
    'Units','pixels', ...
    'Position',[0 0 800 500], ...
    'Renderer','painters', ...
    'Color','w');


%% ------------------------------------------------------------------------
% Panel A: Brain schematic (positions 1-2)
% -------------------------------------------------------------------------

axBrain = subplot(2,4,[1 2]);

plot_brain_regions(EC,brainsurface,cols);

axis(axBrain,'off');


%% ------------------------------------------------------------------------
% Panel B: Example severity vs model parameter relationship
% -------------------------------------------------------------------------

axScatter = subplot(2,4,3);

rng(1);

x = rand(30,1);
y = 0.5*x + randn(30,1)*0.1;

scatter(axScatter,x,y,35,'filled');

xlabel(axScatter,'Depression Severity');
ylabel(axScatter,'Model Parameter');

box(axScatter,'off');
axis(axScatter,'square');

set(axScatter, ...
    'XTick',[], ...
    'YTick',[], ...
    'FontSize',10);

% Schematic linear regression
p = polyfit(x,y,1);
yfit = polyval(p,x);

[~,indmin] = min(x);
[~,indmax] = max(x);

hold(axScatter,'on');

plot(axScatter, ...
    [x(indmin),x(indmax)], ...
    [yfit(indmin),yfit(indmax)], ...
    'r:', ...
    'LineWidth',2);

hold(axScatter,'off');


%% ------------------------------------------------------------------------
% Panel C: Energy landscape
% -------------------------------------------------------------------------

axLandscape = subplot(2,4,4);

[X,Y,Z] = peaks(50);

surf(axLandscape,X,Y,0.5*Z, ...
    'EdgeColor','k', ...
    'LineWidth',0.25);

colormap(axLandscape,redblue(256));
clim(axLandscape,[-6 6]);

shading(axLandscape,'interp');

view(axLandscape,300,25);
axis(axLandscape,'tight');
grid(axLandscape,'off');

set(axLandscape, ...
    'XTick',[], ...
    'YTick',[], ...
    'ZTick',[], ...
    'XColor','none', ...
    'YColor','none', ...
    'ZColor','none');

alpha(axLandscape,0.6);

% Energy landscape colorbar
c = colorbar(axLandscape);

c.Position = [0.93 0.5550 0.0215 0.3900];
c.Ticks = [-4 4];
c.TickLabels = {'Low Energy','High Energy'};
c.Ruler.TickLabelRotation = 90;


%% ------------------------------------------------------------------------
% Panel D: Band-filtered signal
% -------------------------------------------------------------------------

axFiltered = subplot(2,4,5);

hold(axFiltered,'on');

for i = 1:8

    % Cycle through the four prefrontal region colors
    colorIndex = mod(i-1,4) + 2;

    plot(axFiltered, ...
        datNorm(:,i)+i, ...
        'Color',cols(colorIndex,:), ...
        'LineWidth',1);

end

box(axFiltered,'off');

set(axFiltered, ...
    'YTick',[], ...
    'XTick',[], ...
    'YColor','none', ...
    'FontSize',10);

title(axFiltered,{'Band Filtered','Signal Amplitude'});

hold(axFiltered,'off');


%% ------------------------------------------------------------------------
% Panel E: Binarized signal
% -------------------------------------------------------------------------

axBinary = subplot(2,4,6);

hold(axBinary,'on');

for i = 1:8

    chan = datNorm(:,i);

    % Median binarization
    binchan = 2*(chan > median(chan))-1;
    binchan(chan == median(chan)) = -1;

    % Scale for visualization
    binchan = binchan/4 + 0.5;

    colorIndex = mod(i-1,4) + 2;

    plot(axBinary, ...
        binchan+i, ...
        'Color',cols(colorIndex,:), ...
        'LineWidth',1);

end

box(axBinary,'off');

set(axBinary, ...
    'YTick',[], ...
    'XTick',[], ...
    'YColor','none', ...
    'FontSize',10);

title(axBinary,'Binarized Signal');

hold(axBinary,'off');

% Same limits for both signal examples
linkaxes([axFiltered axBinary],'xy');


%% ------------------------------------------------------------------------
% Panel F: Model R2
% -------------------------------------------------------------------------

axR2 = subplot(2,4,7);

hold(axR2,'on');

violinplot( ...
    bandgroups, ...
    R2pmemAll, ...
    'FaceColor','none', ...
    'EdgeColor','k', ...
    'LineWidth',1);

swarmchart( ...
    bandgroups, ...
    R2pmemAll, ...
    8, ...
    rgbMatrix, ...
    'filled', ...
    'MarkerFaceAlpha',0.7);

% Inclusion threshold
yline(axR2,0.6,'--', ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5);

title(axR2,'R^2 value');
xlabel(axR2,'Bands');
ylabel(axR2,'R^2');

box(axR2,'off');

xticks(axR2,1:5);
xticklabels(axR2,allbands);

ylim(axR2,[0 1]);

hold(axR2,'off');


%% ------------------------------------------------------------------------
% Panel G: Model accuracy (r_D)
% -------------------------------------------------------------------------

axAccuracy = subplot(2,4,8);

hold(axAccuracy,'on');

violinplot( ...
    bandgroups, ...
    accuracyAll, ...
    'FaceColor','none', ...
    'EdgeColor','k', ...
    'LineWidth',1);

swarmchart( ...
    bandgroups, ...
    accuracyAll, ...
    8, ...
    rgbMatrix, ...
    'filled', ...
    'MarkerFaceAlpha',0.7);

title(axAccuracy,'r_D value');
xlabel(axAccuracy,'Bands');
ylabel(axAccuracy,'r_D');

box(axAccuracy,'off');

xticks(axAccuracy,1:5);
xticklabels(axAccuracy,allbands);

ylim(axAccuracy,[0 1]);

% Patient legend
h = gobjects(nPts,1);

for pt = 1:nPts

    h(pt) = scatter(axAccuracy, ...
        NaN,NaN,20,cmapPts(pt,:),'filled');

end

legendLabels = arrayfun(@(i) sprintf('S%d',i), ...
    1:nPts,'UniformOutput',false);

legend(axAccuracy,h,legendLabels, ...
    'Location','southwest', ...
    'NumColumns',2);

legend(axAccuracy,'boxoff');

hold(axAccuracy,'off');


%% ========================================================================
% Figure 2: Final panel positions
%
% Retain the manually specified positions from the original script.
% Removed intermediate position calculations that were overwritten.
% ========================================================================

% Top row
axBrain.Position     = [0.0100 0.555 0.500 0.390];
axScatter.Position   = [0.5050 0.555 0.205 0.390];
axLandscape.Position = [0.7055 0.555 0.215 0.390];

% Bottom row
axFiltered.Position  = [0.055 0.105 0.205 0.390];
axBinary.Position    = [0.275 0.105 0.205 0.390];
axR2.Position        = [0.545 0.105 0.205 0.390];
axAccuracy.Position  = [0.795 0.105 0.205 0.390];

set(fig,'Color','w');

drawnow;


%% ========================================================================
% Figure 2: Export individual panels
%
% PNG: Brain, Energy Landscape
% SVG: All remaining panels
% ========================================================================

% Panel A: Brain
exportgraphics(axBrain, ...
    betterLoadSave(fullfile(outDir,'Fig2_Brain.png')), ...
    'Resolution',600);

% Panel B: Model parameter
exportgraphics(axScatter, ...
    betterLoadSave(fullfile(outDir,'Fig2_ModelParameter.svg')), ...
    'ContentType','vector');

% Panel C: Energy landscape
exportgraphics(axLandscape, ...
    betterLoadSave(fullfile(outDir,'Fig2_EnergyLandscape.png')), ...
    'Resolution',600);

% Panel D: Band-filtered signal
exportgraphics(axFiltered, ...
    betterLoadSave(fullfile(outDir,'Fig2_BandFiltered.svg')), ...
    'ContentType','vector');

% Panel E: Binarized signal
exportgraphics(axBinary, ...
    betterLoadSave(fullfile(outDir,'Fig2_Binarized.svg')), ...
    'ContentType','vector');

% Panel F: R2
exportgraphics(axR2, ...
    betterLoadSave(fullfile(outDir,'Fig2_R2.svg')), ...
    'ContentType','vector');

% Panel G: Accuracy
exportgraphics(axAccuracy, ...
    betterLoadSave(fullfile(outDir,'Fig2_Accuracy.svg')), ...
    'ContentType','vector');

% Leave assembled figure open for inspection.
% close(fig);

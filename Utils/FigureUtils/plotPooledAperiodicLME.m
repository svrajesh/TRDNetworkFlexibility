function lme = plotPooledAperiodicLME( ...
    parentAx, ...
    aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,reggroups,recordinggroups, ...
    R2pmemAll,allpts,allregions,alllabels,panelMetric)

%% Settings
ptSize = 45;
lnWidth = 1;
fontSize = 9;
titleFontSize = 12;

R2thr = 0.6;
betaBand = 4;

nPatients = numel(allpts);
cmap = lines(nPatients);

%% =========================================================
% Select metric
% =========================================================
switch panelMetric

    case 'Temperature'
        Y = tempsAll;

    case 'Dwell Time'
        Y = dtAll;

    case 'Transition Probability'
        Y = transPall;

    case 'Entropy Rate'
        Y = entAll;

    otherwise
        error('Invalid panelMetric');

end

%% =========================================================
% Create layout
% =========================================================
tl = tiledlayout( ...
    parentAx,1,1, ...
    'Padding','tight', ...
    'TileSpacing','tight');

%% =========================================================
% Compute pooled X
%
% Beta band only.
%
% For each recording:
%   X = mean(left vmPFC aperiodic slope,
%            right vmPFC aperiodic slope)
%
%   Y = mean(left vmPFC metric,
%            right vmPFC metric)
%
% Only recordings where BOTH vmPFC regions satisfy the
% R2 threshold and have negative aperiodic slopes are used.
% =========================================================
Xpool = [];
Ypool = [];
pid   = [];
rid   = [];

for p = 1:nPatients

    Lpass = find( ...
        ptgroups == p & ...
        bandgroups == betaBand & ...
        reggroups == 1 & ...
        R2pmemAll >= R2thr & ...
        aperiodicsAll < 0);

    Rpass = find( ...
        ptgroups == p & ...
        bandgroups == betaBand & ...
        reggroups == 5 & ...
        R2pmemAll >= R2thr & ...
        aperiodicsAll < 0);

    % Find recordings present in both hemispheres
    [commonRec,idxL,idxR] = intersect( ...
        recordinggroups(Lpass), ...
        recordinggroups(Rpass));

    if isempty(commonRec)
        continue
    end

    Lidx = Lpass(idxL);
    Ridx = Rpass(idxR);

    % Average left/right vmPFC within recording
    Xthis = 0.5 .* ...
        (aperiodicsAll(Lidx) + aperiodicsAll(Ridx));

    Ythis = 0.5 .* ...
        (Y(Lidx) + Y(Ridx));

    Xpool = [Xpool;Xthis(:)];
    Ypool = [Ypool;Ythis(:)];

    pid = [ ...
        pid; ...
        p*ones(numel(Xthis),1)];

    rid = [ ...
        rid; ...
        commonRec(:)];

end

%% =========================================================
% Remove incomplete observations
% =========================================================
valid = ...
    isfinite(Xpool) & ...
    isfinite(Ypool);

Xraw = Xpool(valid);
Yraw = Ypool(valid);
pid  = pid(valid);
rid  = rid(valid);

%% =========================================================
% Within-patient z-score
%
% Used ONLY for visualization.
% =========================================================
Xz = nan(size(Xraw));
Yz = nan(size(Yraw));

for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) >= 2

        Xz(idx) = zscore(Xraw(idx));
        Yz(idx) = zscore(Yraw(idx));

    end

end

%% Remove observations for which z-scoring failed
valid = ...
    isfinite(Xz) & ...
    isfinite(Yz);

Xz   = Xz(valid);
Yz   = Yz(valid);
Xraw = Xraw(valid);
Yraw = Yraw(valid);
pid  = pid(valid);
rid  = rid(valid);

%% =========================================================
% LME
%
% Fit on RAW data.
%
% Each observation represents one recording after
% left/right vmPFC averaging.
%
% Y ~ X + (1|Patient)
% =========================================================
tbl = table( ...
    Xraw, ...
    Yraw, ...
    categorical(pid), ...
    'VariableNames', ...
    {'X','Y','Patient'});

lme = fitlme( ...
    tbl, ...
    'Y ~ X + (1|Patient)');

%% =========================================================
% Pooled normalized regression
%
% Visualization only.
% =========================================================
mdlPlot = fitlm(Xz,Yz);

xfit = linspace( ...
    min(-3), ...
    max(3), ...
    100)';

[yfit,yCI] = predict( ...
    mdlPlot, ...
    xfit);

rPlot = corr( ...
    Xz,Yz, ...
    'Rows','complete');

%% =========================================================
% Standardized LME beta
% =========================================================
XstdWithin = nan(nPatients,1);
YstdWithin = nan(nPatients,1);
nWithin    = nan(nPatients,1);

for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) >= 2

        XstdWithin(p) = std(Xraw(idx));
        YstdWithin(p) = std(Yraw(idx));
        nWithin(p)    = nnz(idx);

    end

end

validSD = ...
    isfinite(XstdWithin) & ...
    isfinite(YstdWithin) & ...
    nWithin > 1;

SDxWithin = sqrt( ...
    sum( ...
        (nWithin(validSD)-1) .* ...
        XstdWithin(validSD).^2) ./ ...
    sum(nWithin(validSD)-1));

SDyWithin = sqrt( ...
    sum( ...
        (nWithin(validSD)-1) .* ...
        YstdWithin(validSD).^2) ./ ...
    sum(nWithin(validSD)-1));

%% Extract LME coefficient
idxX = strcmp( ...
    lme.CoefficientNames, ...
    'X');

betaRaw = lme.Coefficients.Estimate(idxX);

betaStd = ...
    betaRaw .* SDxWithin ./ SDyWithin;

%% =========================================================
% LME p-value
% =========================================================
pLME = lme.Coefficients.pValue(idxX);

if pLME < 1e-16

    pText = 'p < 10^{-16}';

elseif pLME < 0.01

    tp = ceil(log10(pLME));

    pText = sprintf( ...
        'p < 10^{%d}', ...
        tp);

else

    pText = sprintf( ...
        'p = %.3f', ...
        pLME);

end

%% =========================================================
% Plot
% =========================================================
ax = nexttile(tl);
hold(ax,'on');

%% Individual patient data and regression lines
for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) < 2
        continue
    end

    scatter( ...
        ax, ...
        Xz(idx), ...
        Yz(idx), ...
        ptSize, ...
        cmap(p,:), ...
        'filled', 'MarkerFaceAlpha',0.5);

    mdlPt = fitlm( ...
        Xz(idx), ...
        Yz(idx));

    xfitPt = linspace( ...
        min(Xz(idx)), ...
        max(Xz(idx)), ...
        100)';

    yfitPt = predict( ...
        mdlPt, ...
        xfitPt);

    plot( ...
        ax, ...
        xfitPt, ...
        yfitPt, ...
        '--', ...
        'Color',cmap(p,:), ...
        'LineWidth',lnWidth);

end

%% Pooled 95% CI
fill( ...
    ax, ...
    [xfit;flipud(xfit)], ...
    [yCI(:,1);flipud(yCI(:,2))], ...
    [0.8 0.8 0.8], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.3);

%% Pooled regression
plot( ...
    ax, ...
    xfit, ...
    yfit, ...
    'k', ...
    'LineWidth',lnWidth);

xlim(ax,[-3 3]);

%% =========================================================
% Formatting
% =========================================================
xlabel( ...
    ax, ...
    {'vmPFC Aperiodic Slope','(Z Score)'}, ...
    'FontSize',fontSize);

ylabel( ...
    ax, ...
    {panelMetric,'(Z Score)'}, ...
    'FontSize',fontSize);

%% Statistics text
text( ...
    ax, ...
    0.98,0.98, ...
    { ...
    sprintf('Standardized LME β = %.3f',betaStd), ...
    sprintf('Pooled Pearson r = %.3f',rPlot), ...
    sprintf('LME %s',pText)}, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);

axis(ax,'square');
box(ax,'off');

%% Legend only for dwell-time panel
legstr = { ...
    alllabels{1},'', ...
    alllabels{2},'', ...
    alllabels{3},'', ...
    alllabels{4},'', ...
    alllabels{5},'', ...
    alllabels{6},''};

if strcmp(panelMetric,'Dwell Time')

    legend( ...
        ax, ...
        legstr, ...
        'Location','northwestoutside','Box','off');

end

end
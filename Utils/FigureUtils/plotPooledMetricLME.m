function lme = plotPooledMetricLME( ...
    parentAx, ...
    scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,R2pmemAll,allpts,allbands,alllabels,panelMetric)

%% Settings
ptSize = 45;
lnWidth = 1;
fontSize = 9;
titleFontSize = 12;

R2thr = 0.6;
betaBand = 4;

nPatients = numel(allpts);
cmap = lines(nPatients);

tl = tiledlayout( ...
    parentAx,1,1, ...
    'Padding','tight', ...
    'TileSpacing','tight');

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
% Select beta-band data
% =========================================================
idxInclude = ...
    R2pmemAll >= R2thr & ...
    bandgroups == betaBand;

Xraw = scoresAll(idxInclude);
Yraw = Y(idxInclude);
pid  = ptgroups(idxInclude);

%% =========================================================
% Remove incomplete observations
% =========================================================
valid = ...
    isfinite(Xraw) & ...
    isfinite(Yraw);

Xraw = Xraw(valid);
Yraw = Yraw(valid);
pid  = pid(valid);

%% =========================================================
% Require >2 observations per patient
%
% This matches the main analysis.
% =========================================================
keep = false(size(pid));

for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) > 2
        keep(idx) = true;
    end

end

Xraw = Xraw(keep);
Yraw = Yraw(keep);
pid  = pid(keep);

%% =========================================================
% Center BOTH variables within patient
%
% Each patient's mean is calculated from the recordings
% retained above.
%
% Xwithin:
%   severity relative to patient's own mean severity
%
% Ywithin:
%   metric relative to patient's own mean metric
% =========================================================
Xwithin = nan(size(Xraw));
Ywithin = nan(size(Yraw));

for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) > 2

        Xwithin(idx) = ...
            Xraw(idx) - mean(Xraw(idx));

        Ywithin(idx) = ...
            Yraw(idx) - mean(Yraw(idx));

    end

end

%% =========================================================
% Within-patient z-score
% Visualization only
% =========================================================
Xz = nan(size(Xraw));
Yz = nan(size(Yraw));

for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) > 2

        Xz(idx) = zscore(Xraw(idx));
        Yz(idx) = zscore(Yraw(idx));

    end

end

%% =========================================================
% Remove any remaining invalid observations
% =========================================================
valid = ...
    isfinite(Xwithin) & ...
    isfinite(Ywithin) & ...
    isfinite(Xz) & ...
    isfinite(Yz);

Xwithin = Xwithin(valid);
Ywithin = Ywithin(valid);

Xz      = Xz(valid);
Yz      = Yz(valid);

Xraw    = Xraw(valid);
Yraw    = Yraw(valid);

pid     = pid(valid);

%% =========================================================
% LME
%
% Explicit within-patient model:
%
% Ywithin ~ Xwithin + (1|Patient)
% =========================================================
tbl = table( ...
    Xwithin, ...
    Ywithin, ...
    categorical(pid), ...
    'VariableNames', ...
    {'Xwithin','Ywithin','Patient'});

lme = fitlme( ...
    tbl, ...
    'Ywithin ~ Xwithin + (1|Patient)');

%% =========================================================
% Pooled normalized regression
% Visualization only
% =========================================================
mdlPlot = fitlm(Xz,Yz);

xfit = linspace(min(Xz),max(Xz),100)';

[yfit,yCI] = predict( ...
    mdlPlot, ...
    xfit);

rPlot = corr( ...
    Xz, ...
    Yz, ...
    'Rows','complete');

%% =========================================================
% Standardized LME beta
% =========================================================
XstdWithin = nan(nPatients,1);
YstdWithin = nan(nPatients,1);
nWithin    = nan(nPatients,1);

for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) > 2

        XstdWithin(p) = std(Xraw(idx));
        YstdWithin(p) = std(Yraw(idx));
        nWithin(p)    = nnz(idx);

    end

end

validSD = ...
    isfinite(XstdWithin) & ...
    isfinite(YstdWithin) & ...
    nWithin > 2;

if ~any(validSD)

    error( ...
        'Could not calculate within-patient standard deviations.');

end

SDxWithin = sqrt( ...
    sum((nWithin(validSD)-1) .* ...
        XstdWithin(validSD).^2) / ...
    sum(nWithin(validSD)-1));

SDyWithin = sqrt( ...
    sum((nWithin(validSD)-1) .* ...
        YstdWithin(validSD).^2) / ...
    sum(nWithin(validSD)-1));

idxX = strcmp( ...
    lme.CoefficientNames, ...
    'Xwithin');

betaRaw = ...
    lme.Coefficients.Estimate(idxX);

betaStd = ...
    betaRaw * SDxWithin / SDyWithin;

%% =========================================================
% LME p-value
% =========================================================
pLME = ...
    lme.Coefficients.pValue(idxX);

if pLME < 1e-16

    pText = 'p < 10^{-16}';

elseif pLME < 0.01

    tp = ceil(log10(pLME));

    pText = sprintf( ...
        'p < 10^{%d}',tp);

else

    pText = sprintf( ...
        'p = %.3f',pLME);

end

%% =========================================================
% Plot
% =========================================================
ax = nexttile(tl);

hold(ax,'on');

% Individual patient data and regression lines
for p = 1:nPatients

    idx = pid == p;

    if nnz(idx) <= 2
        continue
    end

    scatter( ...
        ax, ...
        Xz(idx), ...
        Yz(idx), ...
        ptSize, ...
        cmap(p,:), ...
        'filled', ...
        'MarkerFaceAlpha',0.5);

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

% Pooled 95% CI
fill( ...
    ax, ...
    [xfit;flipud(xfit)], ...
    [yCI(:,1);flipud(yCI(:,2))], ...
    [0.8 0.8 0.8], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.3);

% Pooled regression
plot( ...
    ax, ...
    xfit, ...
    yfit, ...
    'k', ...
    'LineWidth',lnWidth);

xlim([-2.5,2.5]);

%% =========================================================
% Formatting
% =========================================================
xlabel( ...
    ax, ...
    {'Depression Severity','(Z Score)'}, ...
    'FontSize',fontSize);

ylabel( ...
    ax, ...
    {panelMetric,'(Z Score)'}, ...
    'FontSize',fontSize);

text( ...
    ax, ...
    0.98,0.98, ...
    { ...
    sprintf('Standardized LME \x03b2 = %.3f',betaStd), ...
    sprintf('Pooled Pearson r = %.3f',rPlot), ...
    sprintf('LME %s',pText)}, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);

axis(ax,'square');
box(ax,'off');

legstr = { ...
    alllabels{1},'', ...
    alllabels{2},'', ...
    alllabels{3},'', ...
    alllabels{4},'', ...
    alllabels{5},'', ...
    alllabels{6},''};

if strcmp(panelMetric,'Dwell Time')

    legend( ...
        legstr, ...
        'Location','northwestoutside', ...
        'Box','off');

end

end
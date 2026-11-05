%% Step02 -- Power analysis of weighted stimulation ROI reconstructions
% Computes band summaries and power-vs-score plots.
% FDR correction excludes the broadband 'full' row, as in the original.

if isempty(gcp('nocreate'))
    parpool(32);
end

root = 'USER_DIR';
inputRoot = betterLoadSave(fullfile(root, 'WeightedReconData', task, pt));
outputRoot = betterLoadSave(fullfile(root, 'Analysis', task, pt, 'PowerAnalysis'));
figureRoot = betterLoadSave(fullfile(root, 'Analysis', task, pt, ...
    'Figures', 'PowerAnalysisFigures', 'Weighted'));

if ~isfolder(outputRoot), mkdir(outputRoot); end
if ~isfolder(figureRoot), mkdir(figureRoot); end

chanNames = {'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC','l-dACC','l-MTL', ...
             'r-vmPFC','r-dmPFC','r-vlPFC','r-dlPFC','r-dACC','r-MTL'};
bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};

%% Calculate power and raw spectra
bandpower = [];
rawspec = [];
scoresall = [];
f = [];

for b = 1:numel(bands)
    inputFolder = fullfile(inputRoot, ['clean_', bands{b}]);
    if strcmp(bands{b}, 'raw')
        [rawspec, f, scoresall(b,:)] = speccalc(inputFolder);
    else
        [bandpower(b,:,:), scoresall(b,:)] = bandpowercalc(inputFolder);
    end
end

save(fullfile(outputRoot, 'powerAnalysisWeightedRecon.mat'), ...
    'bandpower', 'scoresall', 'rawspec', 'f', 'chanNames');

%% Power-versus-stim correlations for each band and ROI
nBands = numel(bands) - 1;  % Excludes 'raw', which supplies spectra
nROI = numel(chanNames);
R_values = nan(nBands, nROI);
P_values = nan(nBands, nROI);

for b = 1:nBands
    fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1]);
    tiledlayout(4, 3);
    scores = scoresall(b,:);

    for r = 1:nROI
        nexttile;
        regionPower = squeeze(bandpower(b,:,r));
        bandData = 20 * log10(regionPower); % Preserve original conversion

        [R, P] = corrcoef(bandData, scores);
        R_values(b,r) = R(1,2);
        P_values(b,r) = P(1,2);

        scatter(scores, bandData, 7, [0.066, 0.443, 0.745], 'filled');
        hold on;
        coeffs = polyfit(scores, bandData, 1);
        xFit = linspace(min(scores), max(scores), 100);
        plot(xFit, polyval(coeffs, xFit), 'r', 'LineWidth', 1.5);

        title(chanNames{r});
        xlabel('Scores');
        ylabel('Band Power (dB)');
        text(0.1, -0.35, sprintf('R^2 = %.2f\np = %.3f', ...
            R_values(b,r)^2, P_values(b,r)), ...
            'Units', 'normalized', 'FontSize', 10);
        hold off;
    end

    sgtitle(bands{b});
    saveas(fig, fullfile(figureRoot, sprintf('%s_power_analysis.png', bands{b})));
    close(fig);
end

%% Correlation heatmap with BH-FDR markers
% Preserve the original choice not to FDR-test the 'full' band.
qValues = nan(size(P_values));
qValues(1:end-1,:) = bhAdjustedP(P_values(1:end-1,:));

fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1]);
imagesc(R_values);
colorbar;
set(gca, 'XTick', 1:nROI, 'XTickLabel', chanNames, ...
    'YTick', 1:nBands, 'YTickLabel', bands(1:nBands));
colormap(redblue);
caxis([-1 1]);
title('Heatmap of R Values');
xlabel('Channels');
ylabel('Bands');
hold on;

for b = 1:nBands
    for r = 1:nROI
        if qValues(b,r) < 0.05
            text(r, b, '*', 'Color', 'k', 'FontSize', 12, ...
                'HorizontalAlignment', 'center');
        end
    end
end

saveas(fig, fullfile(figureRoot, 'power_correlation_heatmap.png'));
close(fig);

%% Raw power spectra, colored by score
rawScores = scoresall(end,:);
scoreMin = min(rawScores);
scoreMax = max(rawScores);
if scoreMax > scoreMin
    scoresScaled = (rawScores - scoreMin) ./ (scoreMax - scoreMin);
else
    scoresScaled = zeros(size(rawScores));
end

cmap = turbo(256);
colorIndices = min(256, max(1, round(255 * scoresScaled) + 1));
colorMatrix = cmap(colorIndices,:);

fig = figure('Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1]);
tiledlayout(4, 3);
colormap(cmap);
axall = gobjects(1, nROI);

for r = 1:nROI
    axall(r) = nexttile;
    rawData = squeeze(rawspec(:,:,r));
    hold on;

    for rec = 1:size(rawData,1)
        freqs = f(rec,:);
        keep = freqs <= 75;
        plot(freqs(keep), 10 * log10(rawData(rec,keep)), ...
            'Color', colorMatrix(rec,:), 'LineWidth', 1.5);
    end

    clim([0 1]);
    cb = colorbar;
    cb.Ticks = [0 1];
    cb.TickLabels = {num2str(scoreMin), num2str(scoreMax)};
    title(chanNames{r});
    xlabel('Frequency (Hz)');
    ylabel('Raw Spectrum (dB)');
    hold off;
end

linkaxes(axall, 'x');
sgtitle('Raw Power Spectrum');
saveas(fig, fullfile(figureRoot, 'spectra.png'));
close(fig);

%% Local helper: Benjamini-Hochberg adjusted p-values
function q = bhAdjustedP(p)
    originalSize = size(p);
    p = p(:);
    q = nan(size(p));
    valid = isfinite(p) & p >= 0 & p <= 1;
    pv = p(valid);
    m = numel(pv);

    if m > 0
        [pSorted, order] = sort(pv);
        qSorted = pSorted .* (m ./ (1:m)');
        qSorted = min(1, flipud(cummin(flipud(qSorted))));
        qValid = zeros(m,1);
        qValid(order) = qSorted;
        q(valid) = qValid;
    end

    q = reshape(q, originalSize);
end

%% Step02 -- Weighted reconstruction power analysis
% Requires pt and task in the workspace.
% Requires: bandpowercalc, speccalc, betterLoadSave, redblue
% (The helpers also use get_info and truncate1min.)

rootDir = '';
dataDir = betterLoadSave(fullfile(rootDir, 'WeightedReconData', task, pt));
analysisDir = betterLoadSave(fullfile(rootDir, 'Analysis', task, pt, 'PowerAnalysis'));
figureDir = betterLoadSave(fullfile(rootDir, 'Analysis', task, pt, ...
    'Figures', 'PowerAnalysisFigures', 'Weighted'));

if ~isfolder(analysisDir), mkdir(analysisDir); end
if ~isfolder(figureDir), mkdir(figureDir); end
if isempty(gcp('nocreate')), parpool(32); end

% Region indices match the 12-region reconstruction convention.
% Check MTL versus MTG labels against the saved region file if needed.
chanNames = {'l-vmPFC','l-dmPFC','l-vlPFC','l-dlPFC','l-dACC','l-MTL', ...
             'r-vmPFC','r-dmPFC','r-vlPFC','r-dlPFC','r-dACC','r-MTL'};
bands = {'delta','theta','alpha','beta','gamma','higamma','full'};
nBands = numel(bands);
nRegions = numel(chanNames);

%% Calculate mean band envelopes and raw-signal Welch PSD
bandpower = [];
scoresall = [];

for b = 1:nBands
    folder = fullfile(dataDir, ['clean_', bands{b}]);
    [thisPower, thisScores] = bandpowercalc(folder);
    thisScores = thisScores(:)';

    if b == 1
        nRec = size(thisPower,1);
        assert(size(thisPower,2) == nRegions, 'Expected %d regions.', nRegions);
        bandpower = nan(nBands,nRec,nRegions);
        scoresall = nan(nBands+1,nRec); % final row is raw
    end
    assert(isequal(size(thisPower),[nRec,nRegions]), ...
        'Unexpected recording/region count in band %s.', bands{b});
    assert(numel(thisScores) == nRec, 'Score count mismatch in %s.', bands{b});
    if b > 1
        assert(isequaln(thisScores,scoresall(1,:)), ...
            'Recording order or severity scores differ in band %s.', bands{b});
    end

    bandpower(b,:,:) = thisPower;
    scoresall(b,:) = thisScores;
end

[rawspec,f,rawScores] = speccalc(fullfile(dataDir,'clean_raw'));
rawScores = rawScores(:)';
assert(size(rawspec,1) == nRec && size(rawspec,3) == nRegions, ...
    'Raw spectral dimensions do not match band-envelope dimensions.');
assert(isequaln(rawScores,scoresall(1,:)), ...
    'Raw recording order or severity scores differ from filtered bands.');
scoresall(end,:) = rawScores;

%% Pearson correlations: severity vs. log mean envelope amplitude
R_values = nan(nBands,nRegions);
P_values = nan(nBands,nRegions);

for b = 1:nBands
    fig = figure('Color','w','Units','normalized','Position',[0 0 1 1]);
    tiledlayout(4,3);
    x = scoresall(b,:)';

    for ch = 1:nRegions
        nexttile;
        amp = reshape(bandpower(b,:,ch),[],1);
        y = 20*log10(amp); % mean envelope amplitude in dB; NOT mean squared power
        valid = isfinite(x) & isfinite(y);

        scatter(x(valid),y(valid),7,[0.066 0.443 0.745],'filled');
        hold on;
        if nnz(valid) >= 3 && numel(unique(x(valid))) > 1 && ...
                numel(unique(y(valid))) > 1
            [R,P] = corrcoef(x(valid),y(valid));
            R_values(b,ch) = R(1,2);
            P_values(b,ch) = P(1,2);
            coeff = polyfit(x(valid),y(valid),1);
            xFit = linspace(min(x(valid)),max(x(valid)),100);
            plot(xFit,polyval(coeff,xFit),'r','LineWidth',1.5);
        end
        title(chanNames{ch});
        xlabel('Severity score');
        ylabel('Mean envelope amplitude (dB)');
        text(0.05,0.05,sprintf('R^2 = %.2f\np = %.3g', ...
            R_values(b,ch)^2,P_values(b,ch)), ...
            'Units','normalized','FontSize',9);
        hold off;
    end

    sgtitle(bands{b});
    saveas(fig,fullfile(figureDir,sprintf('%s_power_analysis.png',bands{b})));
    close(fig);
end

%% FDR correction: 6 canonical bands x 12 regions (72 tests)
% Matches the original script: 'full' is plotted but NOT FDR-tested.
% The original called fdr_BH. This local implementation explicitly returns
% Benjamini-Hochberg-adjusted p-values, avoiding ambiguous fdr_BH outputs.
qValues = nan(size(P_values));
qValues(1:nBands-1,:) = bhAdjust(P_values(1:nBands-1,:));

save(fullfile(analysisDir,'powerAnalysisWeightedRecon.mat'), ...
    'bandpower','scoresall','rawspec','f','chanNames', ...
    'R_values','P_values','qValues','bands','-v7.3');

%% Correlation heatmap
fig = figure('Color','w','Units','normalized','Position',[0 0 1 1]);
imagesc(R_values,[-1 1]);
colorbar;
colormap(redblue);
set(gca,'XTick',1:nRegions,'XTickLabel',chanNames, ...
        'YTick',1:nBands,'YTickLabel',bands);
xlabel('Region'); ylabel('Band');
title('Severity correlation with mean envelope amplitude');
hold on;
for b = 1:nBands-1
    for ch = 1:nRegions
        if isfinite(qValues(b,ch)) && qValues(b,ch) < 0.05
            text(ch,b,'*','Color','k','FontSize',12, ...
                'HorizontalAlignment','center');
        end
    end
end
hold off;
saveas(fig,fullfile(figureDir,'power_correlation_heatmap.png'));
close(fig);

%% Raw spectra, colored by severity score
fig = figure('Color','w','Units','normalized','Position',[0 0 1 1]);
tiledlayout(4,3);
cmap = turbo(256);
scores = scoresall(end,:);
finiteScores = scores(isfinite(scores));
assert(~isempty(finiteScores),'No finite raw severity scores.');
lo = min(finiteScores);
hi = max(finiteScores);
if hi == lo
    colorIdx = repmat(128,1,nRec);
    colorLimits = [lo-0.5,hi+0.5];
else
    colorIdx = max(1,min(256,1+round(255*(scores-lo)/(hi-lo))));
    colorLimits = [lo,hi];
end
axAll = gobjects(1,nRegions);

for ch = 1:nRegions
    ax = nexttile;
    axAll(ch) = ax;
    hold(ax,'on');
    for rec = 1:nRec
        if ~isfinite(scores(rec)), continue; end
        freqs = f(rec,:);
        psd = reshape(rawspec(rec,:,ch),1,[]);
        use = isfinite(freqs) & freqs > 0 & freqs <= 75 & ...
              isfinite(psd) & psd > 0;
        plot(ax,freqs(use),10*log10(psd(use)), ...
            'Color',cmap(colorIdx(rec),:),'LineWidth',1.5);
    end
    colormap(ax,cmap);
    clim(ax,colorLimits);
    colorbar(ax);
    title(ax,chanNames{ch});
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'PSD (dB/Hz)');
    hold(ax,'off');
end
linkaxes(axAll,'x');
sgtitle('Raw Welch power spectral density');
saveas(fig,fullfile(figureDir,'spectra.png'));
close(fig);

%% Local BH procedure returning adjusted p-values (q-values)
function q = bhAdjust(p)
    q = nan(size(p));
    mask = isfinite(p) & p >= 0 & p <= 1;
    pv = p(mask);
    if isempty(pv), return; end
    [ps,order] = sort(pv(:));
    m = numel(ps);
    adjusted = ps .* (m ./ (1:m)');
    adjusted = min(1,flipud(cummin(flipud(adjusted))));
    reordered = zeros(m,1);
    reordered(order) = adjusted;
    q(mask) = reordered;
end

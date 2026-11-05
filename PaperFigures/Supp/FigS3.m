%% FIGURE S3 | Beta-band region-ablation sensitivity analysis
% Original raw-data interaction LME and common scaling factor
% are retained, without additional random effects or centering.


% Gather Data
allpts = {};
alllabels = {};
alltasks = {'Resting'};
allbands = {'beta'};

regions = {'l-vmPFC','l-vlPFC','l-dlPFC','l-dACC', ...
    'r-vmPFC','r-vlPFC','r-dlPFC','r-dACC'};

thresholds = [50];

ptgroups = [];
reggroups = [];
recordinggroups = [];

accuracyAll = [];
R2pmemAll = [];
tempsAll = [];
scoresAll = [];
transPall = [];
dtAll = [];
entAll = [];

% Load all ablation models
for pt = 1:length(allpts)

    for tk = 1:length(alltasks)

        for b = 1:length(allbands)

            for thresh = 1:length(thresholds)

                for reg = 1:length(regions)

                    paramfile = betterLoadSave( ...
                        [USER_DIR,'/Analysis/Ablation/' ...
                        alltasks{tk},'/',allpts{pt},'/Params/reg_', ...
                        num2str(reg),'_pMEM_features_', ...
                        num2str(thresholds(thresh)),'/', ...
                        allbands{b},'_params.mat']);

                    load(paramfile);

                    % Grouping variables
                    nParams = length(params);

                    ptgrouping  = pt * ones(nParams,1);
                    reggrouping = reg * ones(nParams,1);

                    % Each params entry corresponds to one recording.
                    recordinggrouping = (1:nParams)';

                    % Model metrics
                    mets = [params.metrics];
                    mets = [mets.sharedMetrics];

                    thisAcc    = [mets.d];
                    thisR2pmem  = [mets.R2pmem];

                    thisTemps = [params.parameters];
                    thisTemps = [thisTemps.effectiveTemperature];

                    % Transition statistics
                    Rall = nan(nParams,1);
                    meanDwellAll = nan(nParams,1);
                    HrateAll = nan(nParams,1);
                    NeffAll = nan(nParams,1);

                    for pf = 1:nParams

                        Dtrans = ...
                            params(pf).sharedFeatures.Dynamics.Dtrans;

                        [Rall(pf), ...
                            meanDwellAll(pf), ...
                            HrateAll(pf), ...
                            NeffAll(pf)] = ...
                            summarizeTransitions(Dtrans);

                    end

                    thisR  = Rall;
                    thisDT = meanDwellAll;
                    thisHR = HrateAll;

                    thisScores = [params.catdiscore];

                    % Accumulate
                    ptgroups = [ptgroups; ptgrouping(:)];
                    reggroups = [reggroups; reggrouping(:)];
                    recordinggroups = ...
                        [recordinggroups; recordinggrouping(:)];

                    accuracyAll = [accuracyAll; thisAcc(:)];
                    R2pmemAll   = [R2pmemAll; thisR2pmem(:)];

                    scoresAll = [scoresAll; thisScores(:)];
                    tempsAll  = [tempsAll; thisTemps(:)];

                    transPall = [transPall; thisR(:)];
                    dtAll     = [dtAll; thisDT(:)];
                    entAll    = [entAll; thisHR(:)];

                end
            end
        end
    end
end

% Settings for analysis
R2thr = 0.6;

nPatients = length(allpts);
nRegions  = length(regions);

metricNames = { ...
    'Temperature', ...
    'Transition Probability', ...
    'Entropy Rate', ...
    'Dwell Time'};

metricData = { ...
    tempsAll, ...
    transPall, ...
    entAll, ...
    dtAll};

% Figure
% Top row:
%   Patient x Ablated Region correlation heatmaps
% Bottom row:
%   Region x Severity LME forest plots
figure( ...
    'Renderer','painters', ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 16 7.5]);

tl = tiledlayout(2,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Analyze each model parameter
for m = 1:length(metricNames)

    Y = metricData{m};

    % TOP: Patient x Region correlation
    axH = nexttile(tl,m);
    hold(axH,'on');

    Rmat = nan(nRegions,nPatients);
    Pmat = nan(nRegions,nPatients);

    for pt = 1:nPatients

        for reg = 1:nRegions

            idx = ...
                ptgroups == pt & ...
                reggroups == reg & ...
                R2pmemAll >= R2thr;

            X = scoresAll(idx);
            Ythis = Y(idx);

            valid = ...
                isfinite(X) & ...
                isfinite(Ythis);

            X = X(valid);
            Ythis = Ythis(valid);

            if numel(X) >= 2

                [r,pv] = corr( ...
                    X,Ythis, ...
                    'Type','Pearson');

                Rmat(reg,pt) = r;
                Pmat(reg,pt) = pv;

            end
        end
    end

    % FDR correction across all patient x region tests
    pVec = Pmat(:);
    validP = isfinite(pVec);

    pFDRvec = nan(size(pVec));

    if any(validP)

        pFDRvec(validP) = ...
            fdr_BH(pVec(validP),0.05);

    end

    pfdr = reshape( ...
        pFDRvec, ...
        size(Pmat));

    % pcolor with explicit cell boundaries
    xEdges = 0.5:1:(nPatients+0.5);
    yEdges = 0.5:1:(nRegions+0.5);

    Rplot = [Rmat Rmat(:,end); ...
             Rmat(end,:) Rmat(end,end)];

    pcolor( ...
        axH, ...
        xEdges, ...
        yEdges, ...
        Rplot);

    shading(axH,'flat');

    axH.YDir = 'reverse';
    axH.CLim = [-1 1];

    colormap(axH,redblue);

    % Axes
    xlim(axH,[0.5 nPatients+0.5]);
    ylim(axH,[0.5 nRegions+0.5]);

    set(axH, ...
        'XTick',1:nPatients, ...
        'XTickLabel',alllabels, ...
        'YTick',1:nRegions, ...
        'YTickLabel',regions, ...
        'TickLength',[0 0], ...
        'FontSize',fontSize);

    % Labels
    if m == 1

        xlabel( ...
            axH, ...
            'Patient', ...
            'FontSize',fontSize);

        ylabel( ...
            axH, ...
            'Ablated Region', ...
            'FontSize',fontSize);

    end

    title( ...
        axH, ...
        metricNames{m}, ...
        'FontSize',titleFontSize);

    % Significance markers
    for pt = 1:nPatients

        for reg = 1:nRegions

            if ~isfinite(pfdr(reg,pt))
                continue
            end

            if pfdr(reg,pt) < 0.001

                text( ...
                    axH, ...
                    pt,reg,'***', ...
                    'Color','k', ...
                    'FontSize',10, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle');

            elseif pfdr(reg,pt) < 0.01

                text( ...
                    axH, ...
                    pt,reg,'**', ...
                    'Color','k', ...
                    'FontSize',10, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle');

            elseif pfdr(reg,pt) < 0.05

                text( ...
                    axH, ...
                    pt,reg,'*', ...
                    'Color','k', ...
                    'FontSize',10, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle');

            end
        end
    end

    box(axH,'off');

    % BOTTOM: Region x Severity LME
    axF = nexttile(tl,4+m);
    hold(axF,'on');

    % Select observations
    idxInclude = ...
        R2pmemAll >= R2thr & ...
        isfinite(scoresAll) & ...
        isfinite(Y) & ...
        isfinite(ptgroups) & ...
        isfinite(reggroups) & ...
        isfinite(recordinggroups);

    X = scoresAll(idxInclude);
    Ylme = Y(idxInclude);

    patient   = ptgroups(idxInclude);
    region    = reggroups(idxInclude);
    recording = recordinggroups(idxInclude);

    % Fit the specified raw-data LME. Scaling is display-only.
    Xraw = X;
    Yraw = Ylme;

    % Categorical variables
    PatientCat   = categorical(patient);
    RegionCat    = categorical(region);
    RecordingCat = categorical(recording);

    tbl = table( ...
        Xraw, ...
        Yraw, ...
        PatientCat, ...
        RegionCat, ...
        RecordingCat, ...
        'VariableNames', ...
        {'Severity','Y','Patient','Region','Recording'});

    % Region x Severity LME
    % No recording-level random effect is included in this analysis.
    lme = fitlme( ...
        tbl, ...
        'Y ~ Severity*Region + (1|Patient) + (1|Patient:Region)' );

    % Extract fixed effects
    coefNames = lme.CoefficientNames;
    beta = lme.Coefficients.Estimate;
    CovBeta = lme.CoefficientCovariance;

    idxSeverity = ...
        strcmp(coefNames,'Severity');

    betaSeverityRaw = beta(idxSeverity);

    % Reference region
    regionLevels = categories(RegionCat);
    refRegion = str2double(regionLevels{1});

    % Find interaction coefficient
    findInteraction = @(r) find( ...
        strcmp(coefNames,sprintf('Severity:Region_%d',r)) | ...
        strcmp(coefNames,sprintf('Region_%d:Severity',r)), ...
        1);

    % Region-specific RAW severity slopes
    slopeRaw   = nan(nRegions,1);
    slopeRawSE = nan(nRegions,1);

    for r = 1:nRegions

        if r == refRegion

            slopeRaw(r) = betaSeverityRaw;

            slopeRawSE(r) = sqrt( ...
                CovBeta(idxSeverity,idxSeverity));

        else

            idxInt = findInteraction(r);

            if ~isempty(idxInt)

                slopeRaw(r) = ...
                    betaSeverityRaw + beta(idxInt);

                slopeRawSE(r) = sqrt( ...
                    CovBeta(idxSeverity,idxSeverity) + ...
                    CovBeta(idxInt,idxInt) + ...
                    2*CovBeta(idxSeverity,idxInt));

            end
        end
    end

    % Standardize slopes using ONE pooled within-patient SD ratio
    % across all ablation conditions. The common scale preserves
    % the ratios of regional raw slopes; it is not region-specific.
    % betaStd = betaRaw * SD(X)within / SD(Y)within
    % The p-values remain those from the RAW LME.
    XstdWithin = nan(nPatients,1);
    YstdWithin = nan(nPatients,1);
    nWithin    = nan(nPatients,1);

    for p = 1:nPatients

        idx = patient == p;

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

    if any(validSD) && ...
            sum((nWithin(validSD)-1)) > 0 && ...
            all(isfinite([ ...
                XstdWithin(validSD); ...
                YstdWithin(validSD)]))

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

        slope = ...
            slopeRaw .* SDxWithin ./ SDyWithin;

        slopeSE = ...
            slopeRawSE .* SDxWithin ./ SDyWithin;

    else

        slope = slopeRaw;
        slopeSE = slopeRawSE;

    end

    % 95% confidence intervals
    % Scale the raw LME confidence intervals by the same
    % SD ratio.
    tCrit = tinv(0.975,lme.DFE);

    CIraw = [ ...
        slopeRaw - tCrit*slopeRawSE, ...
        slopeRaw + tCrit*slopeRawSE];

    if any(validSD) && ...
            isfinite(SDxWithin) && ...
            isfinite(SDyWithin) && ...
            SDyWithin > 0

        scaleFactor = SDxWithin ./ SDyWithin;

        CI = CIraw .* scaleFactor;

    else

        CI = CIraw;

    end

    % Omnibus Region x Severity interaction
    interactionIdx = [];

    for r = 1:nRegions

        if r == refRegion
            continue
        end

        idxInt = findInteraction(r);

        if ~isempty(idxInt)
            interactionIdx(end+1) = idxInt;
        end
    end

    if isempty(interactionIdx)

        pInteraction = NaN;

    else

        L = zeros( ...
            numel(interactionIdx), ...
            numel(coefNames));

        for k = 1:numel(interactionIdx)
            L(k,interactionIdx(k)) = 1;
        end

        pInteraction = coefTest(lme,L);

    end

    % Forest plot
    y = 1:nRegions;

    xline( ...
        axF, ...
        0, ...
        'k-', ...
        'LineWidth',0.75);

    for r = 1:nRegions

        if ~isfinite(slope(r))
            continue
        end

        plot( ...
            axF, ...
            CI(r,:), ...
            [y(r) y(r)], ...
            '-', ...
            'Color',[0.3 0.3 0.3], ...
            'LineWidth',lnWidth);

        plot( ...
            axF, ...
            slope(r), ...
            y(r), ...
            'ko', ...
            'MarkerFaceColor','k', ...
            'MarkerSize',5);

    end

    set(axF, ...
        'YTick',y, ...
        'YTickLabel',regions, ...
        'YDir','reverse', ...
        'FontSize',fontSize, ...
        'TickLength',[0 0]);

    xlabel( ...
        axF, ...
        'Standardized Severity Slope', ...
        'FontSize',fontSize);

    if m == 1

        ylabel( ...
            axF, ...
            'Ablated Region', ...
            'FontSize',fontSize);

    end

    title( ...
        axF, ...
        'Region × Severity', ...
        'FontSize',titleFontSize);

    ylim(axF,[0.5 nRegions+0.5]);

    box(axF,'off');

    % Omnibus interaction p-value
    if ~isfinite(pInteraction)

        interactionText = 'p = NA';

    elseif pInteraction < 0.001

        interactionText = ...
            sprintf('p = %.2g',pInteraction);

    else

        interactionText = ...
            sprintf('p = %.3f',pInteraction);

    end

    text( ...
        axF, ...
        0.98,0.97, ...
        sprintf('Region × Severity: %s',interactionText), ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',fontSize);

end

% Export
statsDir = [USER_DIR,'/StatTables'];
S3stats = exportFigureS3Statistics( ...
    scoresAll, tempsAll, dtAll, transPall, entAll, ...
    ptgroups, reggroups, recordinggroups, R2pmemAll, ...
    allpts, regions, alllabels, ...
    fullfile(statsDir, 'FigureS3_Statistics.xlsx'));
 exportgraphics(gcf,betterLoadSave([USER_DIR,'/PaperFigures/FigS3_HeatmapsAblation.svg']), 'Resolution',600);

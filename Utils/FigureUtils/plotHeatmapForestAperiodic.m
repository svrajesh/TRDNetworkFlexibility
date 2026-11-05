function lmeAll = plotHeatmapForestAperiodic( ...
    parentAx, ...
    aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,reggroups,recordinggroups,R2pmemAll, ...
    allpts,allregions,alllabels,panelMetric)

%% =========================================================
% Settings
% =========================================================
fontSize = 12;
titleFontSize = 15;

R2thr = 0.6;
betaBand = 4;

nPatients = numel(allpts);
nRegions  = numel(allregions);

% Store one LME per region
lmeAll = cell(nRegions,1);

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
% Create 1 x 2 layout
% =========================================================
tl = tiledlayout( ...
    parentAx, ...
    4,1, ...
    'Padding','tight', ...
    'TileSpacing','tight');

%% =========================================================
% PANEL A: Patient x Region correlation heatmap
%
% X = aperiodic slope
% Y = selected model metric
%
% Only beta-band observations are used.
%% =========================================================

axH = nexttile(tl,1,[2 1]);
hold(axH,'on');

Rmat = nan(nRegions,nPatients);
Pmat = nan(nRegions,nPatients);

for p = 1:nPatients

    for r = 1:nRegions

        idx = ...
            ptgroups == p & ...
            bandgroups == betaBand & ...
            reggroups == r & ...
            R2pmemAll >= R2thr;

        X = aperiodicsAll(idx);
        Ythis = Y(idx);

        valid = ...
            isfinite(X) & ...
            isfinite(Ythis);

        X = X(valid);
        Ythis = Ythis(valid);

        if numel(X) >= 2

            [rval,pval] = corr( ...
                X, ...
                Ythis, ...
                'Type','Pearson');

            Rmat(r,p) = rval;
            Pmat(r,p) = pval;

        end
    end
end

%% =========================================================
% FDR correction across all patient x region tests
% =========================================================

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

%% =========================================================
% Heatmap
% =========================================================

pcolor( ...
    axH, ...
    [Rmat Rmat(:,end); ...
     Rmat(end,:) Rmat(end,end)]);

shading(axH,'flat');
axis(axH,'ij');

colormap(axH,redblue);
axH.CLim = [-1 1];

set(axH, ...
    'XTick',0.5+(1:nPatients), ...
    'XTickLabel',1:nPatients, ...
    'YTick',0.5+(1:nRegions), ...
    'YTickLabel',allregions, ...
    'TickLength',[0 0], ...
    'FontSize',fontSize);
xlim([1,nPatients+1]);
xlabel( ...
    axH, ...
    'Subject ID', ...
    'FontSize',fontSize);

ylabel( ...
    axH, ...
    'Brain Region', ...
    'FontSize',fontSize);



%% =========================================================
% Heatmap significance markers
% =========================================================

for p = 1:nPatients

    for r = 1:nRegions

        if ~isfinite(pfdr(r,p))
            continue
        end

        if pfdr(r,p) < 0.001

            text( ...
                axH, ...
                p+0.5, ...
                r+0.6, ...
                '***', ...
                'HorizontalAlignment','center');

        elseif pfdr(r,p) < 0.01

            text( ...
                axH, ...
                p+0.5, ...
                r+0.6, ...
                '**', ...
                'HorizontalAlignment','center');

        elseif pfdr(r,p) < 0.05

            text( ...
                axH, ...
                p+0.5, ...
                r+0.6, ...
                '*', ...
                'HorizontalAlignment','center');

        end
    end
end

%% =========================================================
% PANEL B: Region-wise LME forest plot
%
% Separate model for each region:
%
%     Y ~ X + (1|Patient)
%
% X = regional aperiodic slope
% Y = selected metric
%
% LME is fit on RAW X and Y.
%
% Reported slopes are standardized using pooled
% within-patient SDs, consistent with the pooled
% aperiodic LME:
%
%     beta_std =
%         beta_raw * SD_X_within / SD_Y_within
%% =========================================================

slopeRaw   = nan(nRegions,1);
slopeRawSE = nan(nRegions,1);
pSlope     = nan(nRegions,1);

slope      = nan(nRegions,1);
slopeSE    = nan(nRegions,1);
CI         = nan(nRegions,2);

for r = 1:nRegions

    %% -----------------------------------------------------
    % Select beta-band observations for this region
    %% -----------------------------------------------------

    idxInclude = ...
        bandgroups == betaBand & ...
        reggroups == r & ...
        R2pmemAll >= R2thr;

    X = aperiodicsAll(idxInclude);
    Yreg = Y(idxInclude);

    patient = ptgroups(idxInclude);
    record  = recordinggroups(idxInclude);

    %% -----------------------------------------------------
    % Remove incomplete observations
    %% -----------------------------------------------------

    valid = ...
        isfinite(X) & ...
        isfinite(Yreg) & ...
        isfinite(patient) & ...
        isfinite(record);

    X       = X(valid);
    Yreg    = Yreg(valid);
    patient = patient(valid);
    record  = record(valid);

    %% -----------------------------------------------------
    % Require sufficient data
    %% -----------------------------------------------------

    if numel(X) < 3 || numel(unique(patient)) < 2
        continue
    end

    %% -----------------------------------------------------
    % Within-patient SDs
    %
    % These are used only to standardize the reported
    % coefficient after fitting the raw-data LME.
    %% -----------------------------------------------------

    XstdWithin = nan(nPatients,1);
    YstdWithin = nan(nPatients,1);
    nWithin    = nan(nPatients,1);

    for p = 1:nPatients

        idxP = patient == p;

        if nnz(idxP) >= 2

            XstdWithin(p) = std(X(idxP));
            YstdWithin(p) = std(Yreg(idxP));
            nWithin(p)    = nnz(idxP);

        end
    end

    validSD = ...
        isfinite(XstdWithin) & ...
        isfinite(YstdWithin) & ...
        nWithin > 1;

    if ~any(validSD)
        continue
    end

    %% -----------------------------------------------------
    % Pooled within-patient SD
    %
    % Same calculation as plotPooledAperiodicLME.
    %% -----------------------------------------------------

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

    %% -----------------------------------------------------
    % Categorical patient variable
    %% -----------------------------------------------------

    PatientCat = categorical(patient);

    %% -----------------------------------------------------
    % Region-specific LME
    %
    % FIT ON RAW DATA.
    %
    % Y ~ X + (1|Patient)
    %% -----------------------------------------------------

    tbl = table( ...
        X, ...
        Yreg, ...
        PatientCat, ...
        'VariableNames', ...
        {'X','Y','Patient'});

    try

        lme = fitlme( ...
            tbl, ...
            'Y ~ X + (1|Patient)');

    catch

        warning( ...
            'LME failed for region %s.', ...
            allregions{r});

        continue

    end

    lmeAll{r} = lme;

    %% -----------------------------------------------------
    % Extract raw slope, SE, and p-value
    %% -----------------------------------------------------

    idxX = strcmp( ...
        lme.CoefficientNames, ...
        'X');

    if ~any(idxX)
        continue
    end

    slopeRaw(r) = ...
        lme.Coefficients.Estimate(idxX);

    slopeRawSE(r) = ...
        lme.Coefficients.SE(idxX);

    pSlope(r) = ...
        lme.Coefficients.pValue(idxX);

    %% -----------------------------------------------------
    % Standardize slope and SE
    %
    % beta_std =
    %     beta_raw * SDxWithin / SDyWithin
    %% -----------------------------------------------------

    standardizationFactor = ...
        SDxWithin / SDyWithin;

    slope(r) = ...
        slopeRaw(r) * standardizationFactor;

    slopeSE(r) = ...
        slopeRawSE(r) * standardizationFactor;

    %% -----------------------------------------------------
    % 95% CI for standardized slope
    %% -----------------------------------------------------

    tCrit = tinv( ...
        0.975, ...
        lme.DFE);

    CI(r,:) = [ ...
        slope(r) - tCrit*slopeSE(r), ...
        slope(r) + tCrit*slopeSE(r)];

end

%% =========================================================
% FDR correction across region-wise LME slopes
%
% Correction is applied to the raw LME p-values across
% all tested regions.
%% =========================================================

pFDR = nan(size(pSlope));

validP = isfinite(pSlope);

if any(validP)

    pFDR(validP) = ...
        fdr_BH(pSlope(validP),0.05);

end

%% =========================================================
% Forest plot
%% =========================================================

axF = nexttile(tl,3);
hold(axF,'on');

y = (1:nRegions)';

%% =========================================================
% Zero reference
% =========================================================

xline( ...
    axF, ...
    0, ...
    '--r', ...
    'LineWidth',1);

%% =========================================================
% Plot each regional estimate
%% =========================================================

for r = 1:nRegions

    if ~isfinite(slope(r))
        continue
    end

    % Significant after FDR
    isSig = ...
        isfinite(pFDR(r)) && ...
        pFDR(r) < 0.05;

    if isSig

        lineColor = [0 0 0];

    else

        lineColor = [0.5 0.5 0.5];

    end

    %% Confidence interval

    plot( ...
        axF, ...
        CI(r,:), ...
        [y(r) y(r)], ...
        '-', ...
        'Color',lineColor, ...
        'LineWidth',1.5);

    %% Point estimate

    plot( ...
        axF, ...
        slope(r), ...
        y(r), ...
        'o', ...
        'Color',lineColor, ...
        'MarkerFaceColor',lineColor, ...
        'MarkerSize',5);

end

%% =========================================================
% Formatting
% =========================================================

set(axF, ...
    'YTick',y, ...
    'YTickLabel',allregions, ...
    'YDir','reverse', ...
    'FontSize',fontSize);

xlabel( ...
    axF, ...
    'Standardized LME Slope', ...
    'FontSize',fontSize);

ylabel( ...
    axF, ...
    'Brain Region', ...
    'FontSize',fontSize);


ylim(axF,[0.5 nRegions+0.5]);

box(axF,'off');

%% =========================================================
% Omnibus summary text
%
% There is no single omnibus Region × X interaction here,
% because separate models are fit for each region.
%
% Report the number of regions with significant
% FDR-corrected slopes.
%% =========================================================

nSig = sum( ...
    isfinite(pFDR) & ...
    pFDR < 0.05);

text( ...
    axF, ...
    0.98, ...
    0.98, ...
    sprintf( ...
        '%d/%d regions significant (FDR < 0.05)', ...
        nSig, ...
        nRegions), ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);

end
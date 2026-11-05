function lme = plotBandsHeatmapForest( ...
    parentAx, ...
    scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll, ...
    allpts,allbands,alllabels,panelMetric)

%% Settings
fontSize = 12;
titleFontSize = 15;

R2thr = 0.6;
nPatients = numel(allpts);
nBands    = numel(allbands);

% Beta band and planned comparisons
betaBand = 4;
betaCompareBands = [1 2 3 5]';

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
    parentAx,4,1, ...
    'Padding','tight', ...
    'TileSpacing','tight');

%% =========================================================
% PANEL A: Patient x Band correlation heatmap
% =========================================================
axH = nexttile(tl,1,[2 1]);
hold(axH,'on');

Rmat = nan(nBands,nPatients);
Pmat = nan(nBands,nPatients);

for p = 1:nPatients

    for b = 1:nBands

        idx = ...
            ptgroups == p & ...
            bandgroups == b & ...
            R2pmemAll >= R2thr;

        Xthis = scoresAll(idx);
        Ythis = Y(idx);

        valid = isfinite(Xthis) & isfinite(Ythis);

        Xthis = Xthis(valid);
        Ythis = Ythis(valid);

        if numel(Xthis) > 2

            [r,pv] = corr( ...
                Xthis,Ythis, ...
                'Type','Pearson');

            Rmat(b,p) = r;
            Pmat(b,p) = pv;

        end
    end
end

%% FDR correction
pVec = Pmat(:);
validP = isfinite(pVec);

pFDRvec = nan(size(pVec));

if any(validP)

    pFDRvec(validP) = fdr_BH( ...
        pVec(validP), ...
        0.05);

end

pfdr = reshape( ...
    pFDRvec, ...
    size(Pmat));

%% Heatmap
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
    'YTick',0.5+(1:nBands), ...
    'YTickLabel',{'δ','θ','α','β','γ'}, ...
    'TickLength',[0 0], ...
    'FontSize',fontSize);

xlim([1,nPatients+1]);

xlabel( ...
    axH, ...
    'Subject', ...
    'FontSize',fontSize);

ylabel( ...
    axH, ...
    'Frequency Band', ...
    'FontSize',fontSize);

%% Significance markers
for p = 1:nPatients

    for b = 1:nBands

        if ~isfinite(pfdr(b,p))
            continue
        end

        if pfdr(b,p) < 0.001

            text( ...
                axH, ...
                p+0.5,b+0.6,'***', ...
                'HorizontalAlignment','center');

        elseif pfdr(b,p) < 0.01

            text( ...
                axH, ...
                p+0.5,b+0.6,'**', ...
                'HorizontalAlignment','center');

        elseif pfdr(b,p) < 0.05

            text( ...
                axH, ...
                p+0.5,b+0.6,'*', ...
                'HorizontalAlignment','center');

        end
    end
end

%% =========================================================
% PANEL B: Band x Severity interaction LME
% =========================================================

% ---------------------------------------------------------
% Start with observations surviving R2 threshold
% ---------------------------------------------------------
idxInclude = R2pmemAll >= R2thr;

X = scoresAll(idxInclude);
Ylme = Y(idxInclude);

patient = ptgroups(idxInclude);
band    = bandgroups(idxInclude);
record  = recordinggroups(idxInclude);

% ---------------------------------------------------------
% Remove incomplete observations
% ---------------------------------------------------------
valid = ...
    isfinite(X) & ...
    isfinite(Ylme) & ...
    isfinite(patient) & ...
    isfinite(band) & ...
    isfinite(record);

X       = X(valid);
Ylme    = Ylme(valid);
patient = patient(valid);
band    = band(valid);
record  = record(valid);

%% =========================================================
% Verify severity is constant within each recording
% =========================================================
recordGroups = findgroups(patient,record);

nUniqueX = splitapply( ...
    @(x) numel(unique(x)), ...
    X, ...
    recordGroups);

if any(nUniqueX > 1)

    error( ...
        'Severity is not constant within at least one recording.');

end

%% =========================================================
% Patient x Band inclusion
%
% Exclusion is performed independently for each patient-band
% combination BEFORE centering either variable.
% =========================================================
keep = false(size(patient));

for p = 1:nPatients

    for b = 1:nBands

        idxPB = ...
            patient == p & ...
            band == b;

        % Require >2 valid recordings
        if nnz(idxPB) > 2

            keep(idxPB) = true;

        end
    end
end

if ~any(keep)

    error( ...
        'No observations remain after requiring >2 observations per patient and band.');

end

%% =========================================================
% Apply patient-band-specific exclusions
% =========================================================
Xkeep       = X(keep);
Ykeep       = Ylme(keep);
patientKeep = patient(keep);
bandKeep    = band(keep);
recordKeep  = record(keep);

%% =========================================================
% CENTER BOTH X AND Y WITHIN PATIENT x BAND
%
% Each patient-band combination gets its own reference level.
%
% Xwithin:
%   severity relative to that patient's mean severity
%   across retained recordings in that band.
%
% Ywithin:
%   metric relative to that patient's mean metric
%   across retained recordings in that band.
%
% Thus the LME estimates:
%
%   Does being more severe than usual for this patient
%   correspond to being more/less flexible than usual
%   for this patient?
% =========================================================

patientBandGroups = findgroups( ...
    patientKeep, ...
    bandKeep);

nPBG = max(patientBandGroups);

XpatientBandMean = accumarray( ...
    patientBandGroups, ...
    Xkeep, ...
    [nPBG 1], ...
    @mean, ...
    NaN);

YpatientBandMean = accumarray( ...
    patientBandGroups, ...
    Ykeep, ...
    [nPBG 1], ...
    @mean, ...
    NaN);

Xwithin = ...
    Xkeep - XpatientBandMean(patientBandGroups);

Ywithin = ...
    Ykeep - YpatientBandMean(patientBandGroups);

% Final analysis variables
Ylme    = Ywithin;
patient = patientKeep;
band    = bandKeep;
record  = recordKeep;

%% =========================================================
% Sanity check:
% both variables should have mean ~0 within every patient-band
% =========================================================
checkGroups = findgroups(patient,band);

checkXMeans = splitapply( ...
    @mean, ...
    Xwithin, ...
    checkGroups);

checkYMeans = splitapply( ...
    @mean, ...
    Ylme, ...
    checkGroups);

if any(abs(checkXMeans) > 1e-10)

    error( ...
        'Patient-band centering of X failed.');

end

if any(abs(checkYMeans) > 1e-10)

    error( ...
        'Patient-band centering of Y failed.');

end

%% =========================================================
% Within-patient standard deviations
%
% Used ONLY to standardize final LME slopes.
% =========================================================
XstdWithin = nan(nPatients,1);
YstdWithin = nan(nPatients,1);
nWithin    = nan(nPatients,1);

for p = 1:nPatients

    idx = patient == p;

    if nnz(idx) > 2

        XstdWithin(p) = std(Xwithin(idx));
        YstdWithin(p) = std(Ylme(idx));
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
    sum((nWithin(validSD)-1) .* XstdWithin(validSD).^2) / ...
    sum(nWithin(validSD)-1));

SDyWithin = sqrt( ...
    sum((nWithin(validSD)-1) .* YstdWithin(validSD).^2) / ...
    sum(nWithin(validSD)-1));

%% =========================================================
% Categorical variables
% =========================================================
PatientCat = categorical(patient);
BandCat    = categorical(band);
RecordCat  = categorical(record);

%% =========================================================
% Interaction LME
% =========================================================
tbl = table( ...
    Xwithin, ...
    Ylme, ...
    PatientCat, ...
    BandCat, ...
    RecordCat, ...
    'VariableNames', ...
    {'Xwithin','Y','Patient','Band','Recording'});

lme = fitlme( ...
    tbl, ...
    'Y ~ Xwithin*Band + (1|Patient) + (1|Patient:Recording)');

%% =========================================================
% Extract fixed effects
% =========================================================
coefNames = lme.CoefficientNames;
beta      = lme.Coefficients.Estimate;
CovBeta   = lme.CoefficientCovariance;

idxX = strcmp(coefNames,'Xwithin');

if ~any(idxX)
    error('Could not find fixed-effect coefficient Xwithin.');
end

%% =========================================================
% Reference band
% =========================================================
bandLevels = categories(BandCat);
refBand = str2double(bandLevels{1});

%% =========================================================
% Interaction coefficient helper
% =========================================================
findInteraction = @(b) find( ...
    strcmp(coefNames,sprintf('Xwithin:Band_%d',b)) | ...
    strcmp(coefNames,sprintf('Band_%d:Xwithin',b)), ...
    1);

%% =========================================================
% Band-specific slopes
% =========================================================
slopeRaw   = nan(nBands,1);
slopeRawSE = nan(nBands,1);
slopeP     = nan(nBands,1);

for b = 1:nBands

    Lslope = zeros(1,numel(coefNames));

    Lslope(idxX) = 1;

    if b ~= refBand

        idxInt = findInteraction(b);

        if isempty(idxInt)
            continue
        end

        Lslope(idxInt) = 1;

    end

    slopeRaw(b) = Lslope * beta;

    slopeRawSE(b) = sqrt( ...
        Lslope * CovBeta * Lslope');

    slopeP(b) = coefTest( ...
        lme, ...
        Lslope);

end

%% =========================================================
% FDR correction across band-specific slopes
% =========================================================
slopePFDR = nan(size(slopeP));

validSlopeP = isfinite(slopeP);

if any(validSlopeP)

    slopePFDR(validSlopeP) = ...
        fdr_BH(slopeP(validSlopeP),0.05);

end

%% =========================================================
% Standardize slopes
% =========================================================
standardizationFactor = ...
    SDxWithin / SDyWithin;

slope = ...
    slopeRaw * standardizationFactor;

slopeSE = ...
    slopeRawSE * standardizationFactor;

%% =========================================================
% 95% confidence intervals
% =========================================================
tCrit = tinv(0.975,lme.DFE);

CI = [ ...
    slope - tCrit*slopeSE, ...
    slope + tCrit*slopeSE];

%% =========================================================
% Global Band x Severity interaction
% =========================================================
interactionIdx = [];

for b = 1:nBands

    if b == refBand
        continue
    end

    idxInt = findInteraction(b);

    if ~isempty(idxInt)
        interactionIdx(end+1) = idxInt; %#ok<AGROW>
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

%% =========================================================
% Planned beta magnitude comparisons
% =========================================================
betaCompareBands = betaCompareBands(:);

nContrasts = numel(betaCompareBands);

betaContrastRaw   = nan(nContrasts,1);
betaContrastRawSE = nan(nContrasts,1);
betaContrastP     = nan(nContrasts,1);

betaSlopeRaw = slopeRaw(betaBand);

if isfinite(betaSlopeRaw)

    betaSign = sign(betaSlopeRaw);

    if betaSign == 0
        betaSign = 1;
    end

    idxBetaInt = findInteraction(betaBand);

    for k = 1:nContrasts

        b = betaCompareBands(k);

        if b < 1 || b > nBands || ~isfinite(slopeRaw(b))
            continue
        end

        otherSlopeRaw = slopeRaw(b);

        otherSign = sign(otherSlopeRaw);

        if otherSign == 0
            otherSign = 1;
        end

        Lb = zeros(1,numel(coefNames));

        % Beta slope contribution
        Lb(idxX) = Lb(idxX) + betaSign;

        if betaBand ~= refBand

            if isempty(idxBetaInt)
                continue
            end

            Lb(idxBetaInt) = ...
                Lb(idxBetaInt) + betaSign;

        end

        % Comparison-band slope contribution
        Lb(idxX) = ...
            Lb(idxX) - otherSign;

        if b ~= refBand

            idxOtherInt = findInteraction(b);

            if isempty(idxOtherInt)
                continue
            end

            Lb(idxOtherInt) = ...
                Lb(idxOtherInt) - otherSign;

        end

        betaContrastRaw(k) = ...
            Lb * beta;

        betaContrastRawSE(k) = sqrt( ...
            Lb * CovBeta * Lb');

        betaContrastP(k) = ...
            coefTest(lme,Lb);

    end
end

%% =========================================================
% Standardized beta magnitude differences
% =========================================================
betaContrast = ...
    betaContrastRaw * standardizationFactor;

betaContrastSE = ...
    betaContrastRawSE * standardizationFactor;

%% =========================================================
% FDR across planned beta comparisons
% =========================================================
validContrast = isfinite(betaContrastP);

betaContrastPFDR = nan(size(betaContrastP));

if any(validContrast)

    betaContrastPFDR(validContrast) = ...
        fdr_BH( ...
            betaContrastP(validContrast), ...
            0.05);

end

%% =========================================================
% Significant beta magnitude comparisons
% =========================================================
betaSignificant = ...
    isfinite(betaContrast) & ...
    isfinite(betaContrastPFDR) & ...
    betaContrast > 0 & ...
    betaContrastPFDR < 0.05;

sigBands = ...
    betaCompareBands(betaSignificant);

%% =========================================================
% PANEL B: Forest plot
% =========================================================
axF = nexttile(tl,3);
hold(axF,'on');

y = (1:nBands)';

xline( ...
    axF, ...
    0, ...
    '--r', ...
    'LineWidth',1);

for b = 1:nBands

    if ~isfinite(slope(b))
        continue
    end

    if isfinite(slopePFDR(b)) && slopePFDR(b) < 0.05

        lineColor = [0.3 0.3 0.3];
        markerFace = 'k';

    else

        lineColor = [0.7 0.7 0.7];
        markerFace = [0.7 0.7 0.7];

    end

    plot( ...
        axF, ...
        CI(b,:), ...
        [y(b) y(b)], ...
        '-', ...
        'Color',lineColor, ...
        'LineWidth',1);

    plot( ...
        axF, ...
        slope(b), ...
        y(b), ...
        'ko', ...
        'MarkerFaceColor',markerFace, ...
        'MarkerSize',5);

end

set(axF, ...
    'YTick',y, ...
    'YTickLabel',{'δ','θ','α','β','γ'}, ...
    'YDir','reverse', ...
    'FontSize',fontSize);

xlabel( ...
    axF, ...
    'Standardized Severity Slope', ...
    'FontSize',fontSize);

ylabel( ...
    axF, ...
    'Frequency Band', ...
    'FontSize',fontSize);

ylim([0.5,5.5]);

%% =========================================================
% Global interaction text
% =========================================================
if ~isfinite(pInteraction)

    interactionText = 'p = NA';

elseif pInteraction < 1e-16

    interactionText = 'p < 10^{-16}';

elseif pInteraction < 0.01

    tp = ceil(log10(pInteraction));

    interactionText = sprintf( ...
        'p < 10^{%d}',tp);

else

    interactionText = sprintf( ...
        'p = %.3f',pInteraction);

end

text( ...
    axF, ...
    0.98,0.2, ...
    sprintf('Interaction: %s',interactionText), ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);

%% =========================================================
% Beta magnitude comparison text
% =========================================================
bandSymbols = {'δ','θ','α','β','γ'};

if isempty(sigBands)

    betaText = sprintf( ...
        '|%s| > none', ...
        bandSymbols{betaBand});

else

    sigText = strjoin( ...
        bandSymbols(sigBands), ...
        ', ');

    betaText = sprintf( ...
        '|%s| > %s', ...
        bandSymbols{betaBand}, ...
        sigText);

end

text( ...
    axF, ...
    0.98,0.1, ...
    betaText, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);

box(axF,'off');

end
function results = exportFigureS3Statistics(scoresAll, tempsAll, dtAll, ...
    transPall, entAll, ptgroups, reggroups, recordinggroups, R2pmemAll, ...
    allpts, regions, alllabels, outFile)
% EXPORTFIGURES3STATISTICS  Detailed statistics for Figure S3 (ablation).
%
% Re-fits the FigS3 beta-band raw-severity mixed models EXACTLY as specified:
%   Y ~ Severity*Region + (1|Patient) + (1|Patient:Region)
% No recording random effect, no severity centering and no extra filtering.
% R2pmem >= 0.6 inclusion is applied independently for each metric.
%
% Reports, for each of 4 metrics:
%   * Global Severity x Ablated Region interaction (the figure's p-value).
%   * Eight region-specific severity slopes (raw and commonly standardized),
%     their model-covariance SEs, 95% confidence intervals and Wald p-values.
%   * All original fixed-effect coefficients, with SEs, CIs and p-values.
%   * Participant x Ablated Region Pearson correlations, with BH-FDR
%     separately over the finite correlations for each metric, as in FigS3.
%
% Standardization is ONE shared within-participant SD(X)/SD(Y) ratio across
% all ablated regions for each metric, as implemented in FigS3. The model
% and its p-values are fitted on raw (unstandardized) data.
%
% INPUTS are existing FigS3 workspace variables, in the order above.
% OUTPUT results is written to one FigureS3_Statistics worksheet.
%
% The region-slope tests and fixed-effect rows provide added diagnostics;
% these are not additional significance claims made by the figure.

if nargin < 13 || isempty(outFile)
    outFile = fullfile(pwd, 'FigureS3_Statistics.xlsx');
end

metricNames = {'Temperature','Transition Probability', ...
               'Entropy Rate','Dwell Time'};
metricData = {tempsAll, transPall, entAll, dtAll};
nPatients = numel(allpts);
nRegions = numel(regions);
R2thr = 0.6;

assert(numel(alllabels) == nPatients, ...
    'alllabels and allpts must have the same number of entries.');
assert(nRegions == 8, 'The supplied FigS3 code expects exactly 8 ablations.');

Xall = scoresAll(:);
Pall = ptgroups(:);
Gall = reggroups(:);
RecAll = recordinggroups(:);
Qall = R2pmemAll(:);
for k = 1:4
    assert(numel(metricData{k}) == numel(Xall), ...
        'Metric %s has an unexpected number of observations.', metricNames{k});
end
assert(all([numel(Pall),numel(Gall),numel(RecAll),numel(Qall)] == numel(Xall)), ...
    'Input grouping/accuracy vectors differ in length.');

cols = {'Metric','Analysis','Participant_ID','Ablated_Region', ...
    'Comparison','N_Obs','N_Participants','Estimate','SE','CI_Lower', ...
    'CI_Upper','Statistic','Statistic_Type','DF1','DF2','P_Raw', ...
    'P_FDR','FDR_Family','Model','Notes','Standardized_Beta', ...
    'Standardized_SE','Standardized_CI_Lower','Standardized_CI_Upper'};
rows = cell(0, numel(cols));
modelText = 'Y ~ Severity*Region + (1|Patient) + (1|Patient:Region)';

for m = 1:4
    Yall = metricData{m}(:);
    metric = metricNames{m};

    % ------------------------------------------------------------------
    % 1. Raw-data LME: identical inclusion, variables and formula to FigS3
    % ------------------------------------------------------------------
    idx = Qall >= R2thr & isfinite(Xall) & isfinite(Yall) & ...
        isfinite(Pall) & isfinite(Gall) & isfinite(RecAll);
    x = Xall(idx);
    y = Yall(idx);
    p = Pall(idx);
    g = Gall(idx);
    r = RecAll(idx);
    assert(~isempty(x), 'No LME observations for metric %s.', metric);

    PatientCat   = categorical(p);
    RegionCat    = categorical(g);
    RecordingCat = categorical(r);
    tbl = table(x, y, PatientCat, RegionCat, RecordingCat, ...
        'VariableNames',{'Severity','Y','Patient','Region','Recording'});
    lme = fitlme(tbl, modelText);

    names = lme.CoefficientNames;
    if isstring(names), names = cellstr(names); end
    beta = lme.Coefficients.Estimate;
    V = lme.CoefficientCovariance;
    nC = numel(names);
    nObs = height(tbl);
    nPt = numel(unique(p));
    idxSeverity = find(strcmp(names,'Severity'),1);
    assert(~isempty(idxSeverity), ...
        'Severity coefficient not found for metric %s.', metric);
    levels = categories(RegionCat);
    refRegion = str2double(levels{1});

    % Return coefficient index for each reference-coded interaction term.
    % Matlab may print either variable order in the interaction name.
    interactionIdx = nan(nRegions,1);
    for reg = 1:nRegions
        if reg == refRegion, continue; end
        ii = find(strcmp(names,sprintf('Severity:Region_%d',reg)) | ...
                  strcmp(names,sprintf('Region_%d:Severity',reg)),1);
        assert(~isempty(ii), ...
            'Severity x Region_%d coefficient missing for metric %s.',reg,metric);
        interactionIdx(reg) = ii;
    end

    % Omnibus interaction: same contrast construction as FigS3.
    ixInt = interactionIdx(isfinite(interactionIdx));
    if isempty(ixInt)
        pInteraction = NaN; F = NaN; df1 = NaN; df2 = NaN;
    else
        L = zeros(numel(ixInt),nC);
        for k = 1:numel(ixInt)
            L(k,ixInt(k)) = 1;
        end
        [pInteraction,F,df1,df2] = coefTest(lme,L);
    end
    rows(end+1,:) = makeRow(metric,'Global severity x region interaction', ...
        'All participants','All regions','Joint interaction test',nObs,nPt, ...
        NaN,NaN,NaN,NaN,F,'F',df1,df2,pInteraction,NaN,'None', ...
        modelText,'Joint test of seven region-by-severity interaction terms.'); %#ok<AGROW>

    % ------------------------------------------------------------------
    % 2. Same common within-patient SD ratio used by the Figure S3 plots
    % ------------------------------------------------------------------
    XstdWithin = nan(nPatients,1);
    YstdWithin = nan(nPatients,1);
    nWithin = nan(nPatients,1);
    for ip = 1:nPatients
        ind = p == ip;
        if nnz(ind) >= 2
            XstdWithin(ip) = std(x(ind));
            YstdWithin(ip) = std(y(ind));
            nWithin(ip) = nnz(ind);
        end
    end
    validSD = isfinite(XstdWithin) & isfinite(YstdWithin) & nWithin > 1;
    if any(validSD) && sum(nWithin(validSD)-1) > 0 && ...
            all(isfinite([XstdWithin(validSD);YstdWithin(validSD)]))
        SDxWithin = sqrt(sum((nWithin(validSD)-1) .* ...
            XstdWithin(validSD).^2) / sum(nWithin(validSD)-1));
        SDyWithin = sqrt(sum((nWithin(validSD)-1) .* ...
            YstdWithin(validSD).^2) / sum(nWithin(validSD)-1));
        assert(SDyWithin > 0, ...
            'Pooled within-participant Y SD is zero for metric %s.',metric);
        scale = SDxWithin / SDyWithin;
    else
        % Same fallback as Figure S3 when pooled SD cannot be estimated.
        scale = 1;
    end
    tc = tinv(0.975,lme.DFE);

    % Region-specific slopes. SEs use the FULL fixed-effect covariance,
    % including covariance between Severity and its region interaction.
    for reg = 1:nRegions
        L = zeros(1,nC);
        L(idxSeverity) = 1;
        if reg ~= refRegion
            L(interactionIdx(reg)) = 1;
        end
        est = L*beta;
        se = sqrt(max(0,L*V*L'));
        ciLo = est - tc*se;
        ciHi = est + tc*se;
        [pVal,Fsl,df1sl,df2sl] = coefTest(lme,L);
        signedT = sign(est)*sqrt(Fsl);
        rows(end+1,:) = makeRow(metric,'Ablation-region slope', ...
            'All participants',regions{reg},'Severity slope',nObs,nPt, ...
            est,se,ciLo,ciHi,signedT,'t',df1sl,df2sl,pVal,NaN, ...
            'None (figure does not adjust slope p-values)',modelText, ...
            sprintf('Region %d; common within-patient SD ratio = %.12g; R2 >= %.1f', ...
                reg,scale,R2thr), ...
            est*scale,se*scale,ciLo*scale,ciHi*scale); %#ok<AGROW>
    end

    % Raw fixed effects: additional complete statistics for model audit.
    c = lme.Coefficients;
    for j = 1:height(c)
        tci = tinv(0.975,c.DF(j));
        rows(end+1,:) = makeRow(metric,'Model fixed coefficient', ...
            'All participants','Reference coded',names{j},nObs,nPt, ...
            c.Estimate(j),c.SE(j), ...
            c.Estimate(j)-tci*c.SE(j), ...
            c.Estimate(j)+tci*c.SE(j), ...
            c.tStat(j),'t',1,c.DF(j),c.pValue(j),NaN, ...
            'None',modelText, ...
            'Coefficient shown on original raw-data scale.'); %#ok<AGROW>
    end

    % ------------------------------------------------------------------
    % 3. Heatmap correlations: match FigS3, including n >= 2 and BH family
    % ------------------------------------------------------------------
    Rmat = nan(nRegions,nPatients);
    Pmat = nan(nRegions,nPatients);
    Nmat = zeros(nRegions,nPatients);
    for ip = 1:nPatients
        for reg = 1:nRegions
            ix = Pall == ip & Gall == reg & Qall >= R2thr;
            xc = Xall(ix);
            yc = Yall(ix);
            valid = isfinite(xc) & isfinite(yc);
            xc = xc(valid); yc = yc(valid);
            Nmat(reg,ip) = numel(xc);
            if numel(xc) >= 2
                [rv,pv] = corr(xc,yc,'Type','Pearson');
                Rmat(reg,ip) = rv;
                Pmat(reg,ip) = pv;
            end
        end
    end
    pVec = Pmat(:);
    validP = isfinite(pVec);
    pfdr = nan(size(pVec));
    if any(validP)
        % SAME custom helper used by FigS3, with the same set of tests.
        pfdr(validP) = fdr_BH(pVec(validP),0.05);
    end
    pfdr = reshape(pfdr,size(Pmat));

    for ip = 1:nPatients
        for reg = 1:nRegions
            rv = Rmat(reg,ip);
            pv = Pmat(reg,ip);
            n = Nmat(reg,ip);
            if n < 2, continue; end
            if n > 2 && isfinite(rv) && abs(rv) < 1
                tstat = rv * sqrt((n-2)/(1-rv^2));
            else
                tstat = NaN;
            end
            rows(end+1,:) = makeRow(metric,'Participant-region correlation', ...
                alllabels{ip},regions{reg},'Pearson r',n,1, ...
                rv,NaN,NaN,NaN,tstat,'t',NaN,n-2,pv,pfdr(reg,ip), ...
                sprintf('BH: %s patient x ablation region correlations',metric), ...
                'Pearson corr(Severity,Y)', ...
                'Heatmap test. BH applied to finite p-values within this metric.'); %#ok<AGROW>
        end
    end
end

results = cell2table(rows,'VariableNames',cols);
[outDir,~,~] = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
if exist(outFile,'file'), delete(outFile); end
writetable(results,outFile,'Sheet','FigureS3_Statistics');
fprintf('Figure S3: exported %d rows to %s\n',height(results),outFile);
end

function row = makeRow(varargin)
% 20 core fields plus 4 optional shared-scale slope fields.
row = varargin;
if numel(row) == 20
    row = [row,{NaN,NaN,NaN,NaN}];
elseif numel(row) ~= 24
    error('Expected 20 or 24 reporting fields, received %d.',numel(row));
end
end

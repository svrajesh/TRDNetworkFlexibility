function results = exportFigure3Statistics(scoresAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,recordinggroups,R2pmemAll,allpts,allbands,alllabels,outFile)
% EXPORTFIGURE3STATISTICS  Refit Figure 3 models and export one Excel sheet.
%
% Call after the data-loading portion of Fig3.m:
% results = exportFigure3Statistics(scoresAll,tempsAll,dtAll,transPall, ...
%     entAll,ptgroups,bandgroups,recordinggroups,R2pmemAll, ...
%     allpts,allbands,alllabels,fullfile(outDir,'Figure3_Statistics.xlsx'));
%
% Reproduces the analysis specifications in plotPooledMetricLME and
% plotBandsHeatmapForest. Requires Statistics and Machine Learning Toolbox.
% One worksheet, metric blocks ordered by primary analyses first.
% Participant_ID uses alllabels (e.g. Subject 1), not allpts.
% Confidence intervals are calculated from the contrast covariance and
% test denominator DF rather than calling coefCI (MATLAB compatibility).
%
% IMPORTANT: Beta magnitude comparisons use observed slope signs, as in
% the figure code. Their p-values are conditional on data-chosen signs;
% interpret as exploratory, not prespecified absolute-value tests.

if nargin < 13 || isempty(outFile)
    outFile = fullfile(pwd,'Figure3_Statistics.xlsx');
end
if exist(outFile,'file'), delete(outFile); end

metricNames = {'Temperature','Transition Probability','Entropy Rate','Dwell Time'};
metricValues = {tempsAll,transPall,entAll,dtAll};
R2thr = 0.6;
nPatients = numel(allpts);
if numel(alllabels) ~= nPatients
    error('alllabels must contain one display label per participant.');
end
nBands = numel(allbands);
betaBand = 4;
compareBands = [1 2 3 5];

% A consistent schema allows all tests to occupy a single Excel sheet.
cols = {'Metric','Analysis','Participant_ID','Band','Comparison','N_Obs','N_Participants', ...
    'Estimate','SE','CI_Lower','CI_Upper','Statistic','Statistic_Type', ...
    'DF1','DF2','P_Raw','P_FDR','FDR_Family','Model','Notes', ...
    'Standardized_Beta','Standardized_SE','Standardized_CI_Lower','Standardized_CI_Upper'};
rows = cell(0,numel(cols));

for m = 1:numel(metricNames)
    metric = metricNames{m};
    Yall = metricValues{m}(:);
    Xall = scoresAll(:);
    P = ptgroups(:); B = bandgroups(:); R = recordinggroups(:);
    Q = R2pmemAll(:);
    if ~isequal(numel(Yall),numel(Xall),numel(P),numel(B),numel(R),numel(Q))
        error('Input vector lengths differ for metric %s.',metric);
    end

    %% 1. Beta-only within-patient LME
    idx = Q >= R2thr & B == betaBand & isfinite(Xall) & isfinite(Yall);
    x = Xall(idx); y = Yall(idx); p = P(idx);
    keep = false(size(p));
    for i = 1:nPatients
        keep(p == i) = nnz(p == i) > 2;
    end
    x = x(keep); y = y(keep); p = p(keep);
    if isempty(x), error('No beta observations for %s.',metric); end
    xw = x; yw = y;
    for i = unique(p(:))'
        ix = p == i;
        xw(ix) = x(ix)-mean(x(ix));
        yw(ix) = y(ix)-mean(y(ix));
    end
    T = table(xw,yw,categorical(p),'VariableNames',{'Xwithin','Ywithin','Patient'});
    mdl = fitlme(T,'Ywithin ~ Xwithin + (1|Patient)');
    ix = strcmp(mdl.CoefficientNames,'Xwithin');
    c = mdl.Coefficients;
    % Avoid coefCI: calculate the coefficient-specific t interval directly.
    tcrit = tinv(0.975,c.DF(ix));
    ciLow = c.Estimate(ix)-tcrit*c.SE(ix);
    ciHigh = c.Estimate(ix)+tcrit*c.SE(ix);
    stdFactor = pooledWithinSD(x,p)/pooledWithinSD(y,p);
    rows(end+1,:) = makeRow(metric,'Beta-only LME','All participants',allbands{betaBand},'', ...
        numel(x),numel(unique(p)),c.Estimate(ix),c.SE(ix),ciLow,ciHigh, ...
        c.tStat(ix),'t',c.DF(ix),NaN,c.pValue(ix),NaN,'None', ...
        'Ywithin ~ Xwithin + (1|Patient)', ...
        sprintf('Within-patient centered; standardized beta = %.8g; R2 threshold %.2f', ...
        c.Estimate(ix)*stdFactor,R2thr), ...
        c.Estimate(ix)*stdFactor,c.SE(ix)*stdFactor,ciLow*stdFactor,ciHigh*stdFactor);

    %% 2. Participant x band Pearson correlations; BH within each metric
    corrRows = cell(0,numel(cols));
    for i = 1:nPatients
        for b = 1:nBands
            ix = P == i & B == b & Q >= R2thr & isfinite(Xall) & isfinite(Yall);
            n = nnz(ix);
            if n <= 2, continue; end
            [rv,pv] = corr(Xall(ix),Yall(ix),'Type','Pearson');
            corrRows(end+1,:) = makeRow(metric,'Participant-band correlation',alllabels{i}, ...
                allbands{b},'',n,1,rv,NaN,NaN,NaN, ...
                rv*sqrt((n-2)/max(1-rv^2,eps)),'t',n-2,NaN,pv,NaN, ...
                sprintf('BH: %s participant x band correlations',metric), ...
                'Pearson correlation',sprintf('Complete cases; R2 >= %.2f',R2thr));
        end
    end
    if ~isempty(corrRows)
        q = bhAdjust(cell2mat(corrRows(:,16)));
        corrRows(:,17) = num2cell(q);
        rows = [rows;corrRows]; %#ok<AGROW>
    end

    %% 3. Cross-band within-patient x band interaction LME
    ix = Q >= R2thr & isfinite(Xall) & isfinite(Yall) & ...
        isfinite(P) & isfinite(B) & isfinite(R);
    x = Xall(ix); y = Yall(ix); p = P(ix); b = B(ix); r = R(ix);
    grp = findgroups(p,r);
    nu = splitapply(@(z) numel(unique(z)),x,grp);
    if any(nu > 1)
        error('Severity varies within patient-recording for %s.',metric);
    end
    keep = false(size(p));
    for i = 1:nPatients
        for j = 1:nBands
            ib = p == i & b == j;
            if nnz(ib) > 2, keep(ib) = true; end
        end
    end
    x = x(keep); y = y(keep); p = p(keep); b = b(keep); r = r(keep);
    if isempty(x), error('No cross-band observations for %s.',metric); end
    g = findgroups(p,b);
    % Map each observation to its own patient-band mean.
    xmeans = splitapply(@mean,x,g);
    ymeans = splitapply(@mean,y,g);
    xw = x-xmeans(g);
    yw = y-ymeans(g);
    T = table(xw,yw,categorical(p),categorical(b),categorical(r), ...
        'VariableNames',{'Xwithin','Y','Patient','Band','Recording'});
    mdl = fitlme(T,'Y ~ Xwithin*Band + (1|Patient) + (1|Patient:Recording)');
    coefNames = mdl.CoefficientNames;
    betas = mdl.Coefficients.Estimate;
    V = mdl.CoefficientCovariance;
    idxX = find(strcmp(coefNames,'Xwithin'));
    levels = categories(T.Band);
    refBand = str2double(levels{1});
    scale = pooledWithinSD(xw,p)/pooledWithinSD(yw,p);
    slopeL = cell(nBands,1); slope = nan(nBands,1);
    slopeRows = cell(0,numel(cols));
    slopeP = nan(nBands,1);
    for j = 1:nBands
        L = zeros(1,numel(betas)); L(idxX) = 1;
        if j ~= refBand
            ii = findInteraction(coefNames,j);
            if isempty(ii), continue; end
            L(ii) = 1;
        end
        slopeL{j} = L;
        slope(j) = L*betas;
        se = sqrt(L*V*L');
        [pv,F,df1,df2] = coefTest(mdl,L);
        tstat = sign(L*betas)*sqrt(F);
        tc = tinv(0.975,df2);
        slopeP(j) = pv;
        slopeRows(end+1,:) = makeRow(metric,'Cross-band slope','All participants',allbands{j},'', ...
            numel(x),numel(unique(p)),slope(j),se, ...
            slope(j)-tc*se,slope(j)+tc*se,tstat,'t',df2,NaN,pv,NaN, ...
            sprintf('BH: %s five band slopes',metric), ...
            'Y ~ Xwithin*Band + (1|Patient) + (1|Patient:Recording)', ...
            sprintf('Patient-band centered; standardized slope = %.8g; standardized SE = %.8g', ...
            slope(j)*scale,se*scale), ...
            slope(j)*scale,se*scale,(slope(j)-tc*se)*scale,(slope(j)+tc*se)*scale);
    end
    q = bhAdjust(slopeP);
    for j = 1:size(slopeRows,1)
        bb = find(strcmp(allbands,slopeRows{j,4}),1);
        slopeRows{j,17} = q(bb);
    end
    rows = [rows;slopeRows]; %#ok<AGROW>

    %% 4. Global interaction F test
    ints = [];
    for j = 1:nBands
        if j ~= refBand
            ii = findInteraction(coefNames,j);
            if ~isempty(ii), ints(end+1) = ii; end %#ok<AGROW>
        end
    end
    if ~isempty(ints)
        L = zeros(numel(ints),numel(betas));
        for k = 1:numel(ints), L(k,ints(k)) = 1; end
        [pv,F,df1,df2] = coefTest(mdl,L);
        rows(end+1,:) = makeRow(metric,'Global severity x band interaction','All participants','All','', ...
            numel(x),numel(unique(p)),NaN,NaN,NaN,NaN,F,'F',df1,df2,pv,NaN, ...
            'None','Y ~ Xwithin*Band + (1|Patient) + (1|Patient:Recording)', ...
            'Omnibus joint test of severity-by-band interaction coefficients');
    end

    %% 5. Observed-sign beta magnitude contrasts; BH over four tests
    contrastRows = cell(0,numel(cols));
    for j = compareBands
        if isempty(slopeL{betaBand}) || isempty(slopeL{j}), continue; end
        sb = sign(slope(betaBand)); if sb == 0, sb = 1; end
        sj = sign(slope(j)); if sj == 0, sj = 1; end
        L = sb*slopeL{betaBand}-sj*slopeL{j};
        est = L*betas; se = sqrt(L*V*L');
        [pv,F,~,df2] = coefTest(mdl,L);
        tc = tinv(0.975,df2);
        contrastRows(end+1,:) = makeRow(metric,'Beta magnitude contrast','All participants', ...
            allbands{betaBand},sprintf('|beta| - |%s|',allbands{j}), ...
            numel(x),numel(unique(p)),est,se,est-tc*se,est+tc*se, ...
            sign(est)*sqrt(F),'t',df2,NaN,pv,NaN, ...
            sprintf('BH: %s four beta contrasts',metric), ...
            'Observed-sign linear contrast from cross-band LME', ...
            sprintf('Data-dependent signs; exploratory; standardized difference = %.8g',est*scale), ...
            est*scale,se*scale,(est-tc*se)*scale,(est+tc*se)*scale);
    end
    if ~isempty(contrastRows)
        contrastRows(:,17) = num2cell(bhAdjust(cell2mat(contrastRows(:,16))));
        rows = [rows;contrastRows]; %#ok<AGROW>
    end
end

% Arrange each metric as a continuous block, with primary model results
% first and participant-level correlations last. Keep one rectangular table
% so Excel filters and statistical columns remain aligned.
metricOrder = {'Temperature','Transition Probability','Entropy Rate','Dwell Time'};
analysisOrder = {'Beta-only LME','Global severity x band interaction', ...
    'Cross-band slope','Beta magnitude contrast','Participant-band correlation'};
metricRank = zeros(size(rows,1),1);
analysisRank = zeros(size(rows,1),1);
for k = 1:size(rows,1)
    metricRank(k) = find(strcmp(metricOrder,rows{k,1}),1);
    analysisRank(k) = find(strcmp(analysisOrder,rows{k,2}),1);
end
[~,order] = sortrows([metricRank analysisRank (1:size(rows,1))']);
rows = rows(order,:);
results = cell2table(rows,'VariableNames',cols);
% A single worksheet with one header row and filterable, rectangular data.
writetable(results,outFile,'Sheet','Figure3_Statistics','WriteMode','overwritesheet');
fprintf('Exported %d Figure 3 statistical results to %s\n',height(results),outFile);
end

function row = makeRow(varargin)
% Four optional standardized statistics are appended for LME effects.
% Correlations and omnibus tests retain missing values in these columns.
row = varargin;
if numel(row) == 20
    row = [row, {NaN,NaN,NaN,NaN}];
elseif numel(row) ~= 24
    error('Unexpected reporting row width: %d.',numel(row));
end
end

function idx = findInteraction(names,b)
idx = find(strcmp(names,sprintf('Xwithin:Band_%d',b)) | ...
           strcmp(names,sprintf('Band_%d:Xwithin',b)),1);
end

function sd = pooledWithinSD(y,p)
num = 0; den = 0;
for id = unique(p(:))'
    yy = y(p==id);
    if numel(yy) > 2
        num = num + (numel(yy)-1)*var(yy);
        den = den + numel(yy)-1;
    end
end
sd = sqrt(num/den);
end

function q = bhAdjust(p)
q = nan(size(p));
ix = find(isfinite(p));
if isempty(ix), return; end
[ps,order] = sort(p(ix));
n = numel(ps);
qs = min(1,ps(:).*n./(1:n)');
qs = flipud(cummin(flipud(qs)));
q(ix(order)) = qs;
end

function results = exportFigure5Statistics(aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
    ptgroups,bandgroups,reggroups,recordinggroups,R2pmemAll, ...
    allpts,allregions,alllabels,outFile)
% EXPORTFIGURE5STATISTICS  Recreate the Figure 5 statistical analyses.
%
% Run after the data-loading section in Fig5.m:
%   exportFigure5Statistics(aperiodicsAll,tempsAll,dtAll,transPall,entAll, ...
%       ptgroups,bandgroups,reggroups,recordinggroups,R2pmemAll, ...
%       allpts,allregions,alllabels,fullfile(outDir,'Figure5_Statistics.xlsx'));
%
% Matches plotPooledAperiodicLME and the regional aperiodic heatmap/forest:
%   1. Bilateral vmPFC pooled raw-data LME (matched recordings only).
%   2. Eight regional raw-data LMEs; BH FDR within metric across regions.
%   3. Patient x region Pearson correlations; BH FDR within metric.
% LME estimates are standardized after fitting with pooled within-patient SDs.

if nargin < 14 || isempty(outFile)
    outFile = fullfile(pwd,'Figure5_Statistics.xlsx');
end
if exist(outFile,'file'), delete(outFile); end

metricNames = {'Temperature','Transition Probability','Entropy Rate','Dwell Time'};
metricValues = {tempsAll,transPall,entAll,dtAll};
cols = {'Metric','Analysis','Participant_ID','Band','Comparison','N_Obs','N_Participants', ...
    'Estimate','SE','CI_Lower','CI_Upper','Statistic','Statistic_Type', ...
    'DF1','DF2','P_Raw','P_FDR','FDR_Family','Model','Notes', ...
    'Standardized_Beta','Standardized_SE','Standardized_CI_Lower','Standardized_CI_Upper'};
rows = cell(0,numel(cols));

nPatients = numel(allpts);
nRegions = numel(allregions);
if numel(alllabels) ~= nPatients
    error('One display label is required for every participant.');
end
R2thr = 0.6;
betaBand = 4;
P = ptgroups(:); B = bandgroups(:); R = reggroups(:);
REC = recordinggroups(:); Q = R2pmemAll(:); Xall = aperiodicsAll(:);

for m = 1:numel(metricNames)
    metric = metricNames{m};
    Y = metricValues{m}(:);
    if ~isequal(numel(Xall),numel(Y),numel(P),numel(B),numel(R),numel(REC),numel(Q))
        error('Input vectors have inconsistent lengths.');
    end

    %% Matched left/right vmPFC: reproduce the pooled Figure 5 panel
    Xpool = []; Ypool = []; pid = [];
    for p = 1:nPatients
        Lpass = find(P==p & B==betaBand & R==1 & Q>=R2thr & Xall<0);
        Rpass = find(P==p & B==betaBand & R==5 & Q>=R2thr & Xall<0);
        [~,idxL,idxR] = intersect(REC(Lpass),REC(Rpass));
        if isempty(idxL), continue; end
        left = Lpass(idxL); right = Rpass(idxR);
        Xthis = 0.5*(Xall(left)+Xall(right));
        Ythis = 0.5*(Y(left)+Y(right));
        Xpool = [Xpool;Xthis(:)]; %#ok<AGROW>
        Ypool = [Ypool;Ythis(:)]; %#ok<AGROW>
        pid = [pid;p*ones(numel(Xthis),1)]; %#ok<AGROW>
    end
    keep = isfinite(Xpool) & isfinite(Ypool);
    Xpool = Xpool(keep); Ypool = Ypool(keep); pid = pid(keep);

    % The plotted function also removes patients for whom z-scoring fails.
    % Apply this before fitting the pooled raw-data LME.
    keep = false(size(pid));
    for p = 1:nPatients
        ix = pid==p;
        if nnz(ix)>=2 && std(Xpool(ix))>0 && std(Ypool(ix))>0
            keep(ix)=true;
        end
    end
    Xpool = Xpool(keep); Ypool = Ypool(keep); pid = pid(keep);
    if ~isempty(Xpool)
        tbl = table(Xpool,Ypool,categorical(pid), ...
            'VariableNames',{'X','Y','Patient'});
        mdl = fitlme(tbl,'Y ~ X + (1|Patient)');
        c = mdl.Coefficients;
        ix = strcmp(c.Name,'X');
        scale = withinSD(Xpool,pid,2)/withinSD(Ypool,pid,2);
        rows(end+1,:) = coefficientRow(metric,'Bilateral vmPFC pooled LME', ...
            'All participants','l/r-vmPFC averaged',mdl,ix,numel(Xpool), ...
            numel(unique(pid)),scale,'Matched beta recordings; both slopes < 0 and both R2 >= 0.6; fitted on raw averages.'); %#ok<AGROW>

        Xz = nan(size(Xpool)); Yz = nan(size(Ypool));
        for p = unique(pid(:))'
            ix = pid==p;
            Xz(ix) = zscore(Xpool(ix));
            Yz(ix) = zscore(Ypool(ix));
        end
        [rv,pv] = corr(Xz,Yz,'Rows','complete','Type','Pearson');
        rows(end+1,:) = makeRow(metric,'Pooled vmPFC visual correlation', ...
            'All participants','beta','l/r-vmPFC averaged',numel(Xz), ...
            numel(unique(pid)),rv,NaN,NaN,NaN,NaN,'Pearson r',NaN,NaN,pv,NaN, ...
            'None','Correlation of within-patient z values', ...
            'Descriptive pooled correlation from figure visualization; not the LME significance test.'); %#ok<AGROW>
    end

    %% Regional LMEs (beta band only)
    regionRows = cell(0,numel(cols));
    for r = 1:nRegions
        ix = B==betaBand & R==r & Q>=R2thr & isfinite(Xall) & ...
             isfinite(Y) & isfinite(P) & isfinite(REC);
        x = Xall(ix); y = Y(ix); p = P(ix);
        if numel(x)<3 || numel(unique(p))<2, continue; end
        sdX = withinSD(x,p,2);
        sdY = withinSD(y,p,2);
        if ~isfinite(sdX) || ~isfinite(sdY), continue; end
        tbl = table(x,y,categorical(p),'VariableNames',{'X','Y','Patient'});
        try
            mdl = fitlme(tbl,'Y ~ X + (1|Patient)');
        catch ME
            warning('Figure 5 regional model failed (%s, %s): %s', ...
                metric,allregions{r},ME.message);
            continue;
        end
        idxX = strcmp(mdl.CoefficientNames,'X');
        regionRows(end+1,:) = coefficientRow(metric,'Regional aperiodic LME', ...
            'All participants',allregions{r},mdl,idxX,numel(x), ...
            numel(unique(p)),sdX/sdY, ...
            'Beta band; R2 >= 0.6; raw-data random-intercept LME.'); %#ok<AGROW>
    end
    if ~isempty(regionRows)
        q = bhAdjust(cell2mat(regionRows(:,16)));
        regionRows(:,17) = num2cell(q);
        regionRows(:,18) = repmat({sprintf('BH: %s regional slopes',metric)},size(regionRows,1),1);
        rows = [rows;regionRows]; %#ok<AGROW>
    end

    %% Pearson heatmap: correction across ALL valid patient x region cells
    corrRows = cell(0,numel(cols));
    for p = 1:nPatients
        for r = 1:nRegions
            ix = P==p & B==betaBand & R==r & Q>=R2thr & ...
                 isfinite(Xall) & isfinite(Y);
            n = nnz(ix);
            if n < 2, continue; end
            [rv,pv] = corr(Xall(ix),Y(ix),'Type','Pearson');
            corrRows(end+1,:) = makeRow(metric,'Participant-region correlation', ...
                alllabels{p},'beta',allregions{r},n,1,rv,NaN,NaN,NaN, ...
                NaN,'Pearson r',NaN,NaN,pv,NaN, ...
                sprintf('BH: %s patient x region correlations',metric), ...
                'Pearson correlation','Complete cases; R2 >= 0.6.'); %#ok<AGROW>
        end
    end
    if ~isempty(corrRows)
        corrRows(:,17) = num2cell(bhAdjust(cell2mat(corrRows(:,16))));
        rows = [rows;corrRows]; %#ok<AGROW>
    end
end

% Primary pooled analysis, regional forest results, then heatmap statistics.
orderNames = {'Bilateral vmPFC pooled LME','Regional aperiodic LME', ...
    'Pooled vmPFC visual correlation','Participant-region correlation'};
rank = zeros(size(rows,1),2);
for i = 1:size(rows,1)
    rank(i,1) = find(strcmp(metricNames,rows{i,1}),1);
    rank(i,2) = find(strcmp(orderNames,rows{i,2}),1);
end
[~,idx] = sortrows([rank (1:size(rows,1))']);
results = cell2table(rows(idx,:),'VariableNames',cols);
writetable(results,outFile,'Sheet','Figure5_Statistics','WriteMode','overwritesheet');
fprintf('Exported %d Figure 5 rows to %s\n',height(results),outFile);
end

function row = coefficientRow(metric,analysis,participant,region,mdl,ix,n,nSub,scale,note)
c = mdl.Coefficients;
crit = tinv(0.975,mdl.DFE);
est = c.Estimate(ix); se = c.SE(ix);
lo = est-crit*se; hi = est+crit*se;
row = makeRow(metric,analysis,participant,'beta',region,n,nSub, ...
    est,se,lo,hi,c.tStat(ix),'t',c.DF(ix),NaN,c.pValue(ix),NaN, ...
    'None','Y ~ X + (1|Patient)',note,est*scale,se*scale,lo*scale,hi*scale);
end

function sd = withinSD(x,p,minN)
num = 0; den = 0;
for id = unique(p(:))'
    v = x(p==id);
    if numel(v)>=minN
        num = num + (numel(v)-1)*var(v);
        den = den + numel(v)-1;
    end
end
if den <= 0
    sd = NaN;
else
    sd = sqrt(num/den);
end
end

function q = bhAdjust(p)
q = nan(size(p));
valid = find(isfinite(p));
if isempty(valid), return; end
[ps,order] = sort(p(valid));
n = numel(ps);
qs = min(1,ps(:).*n./(1:n)');
qs = flipud(cummin(flipud(qs)));
q(valid(order)) = qs;
end

function row = makeRow(varargin)
row = varargin;
if numel(row)==20
    row = [row,{NaN,NaN,NaN,NaN}];
elseif numel(row)~=24
    error('Unexpected reporting row width: %d',numel(row));
end
end

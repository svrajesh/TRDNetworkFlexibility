function results = exportFigure4Statistics(T,alllabels,outFile)
% EXPORTFIGURE4STATISTICS  Export the exact models used in Fig4.m.
%
% Run after Fig4.m has constructed the long-form table T:
%   exportFigure4Statistics(T,alllabels,fullfile(outDir,'Figure4_Statistics.xlsx'));
%
% One Excel worksheet containing the five model outcomes, their fixed
% effects, and raw/z-scored ON/OFF descriptive summaries. No FDR correction
% is added because Fig4.m reports unadjusted stimulation p-values.
%
% The displayed ON/OFF values are z-scored within patient. LME coefficients
% are fitted on raw outcomes. A condition coefficient divided by the pooled
% within-patient outcome SD is reported in SD units of the outcome.

if nargin < 3 || isempty(outFile)
    outFile = fullfile(pwd,'Figure4_Statistics.xlsx');
end
if exist(outFile,'file'), delete(outFile); end

required = {'subject','measure','condition','value','accuracy'};
if ~all(ismember(required,T.Properties.VariableNames))
    error('T must contain: %s.',strjoin(required,', '));
end
nPatients = numel(alllabels);
metricNames = {'Temperature','Transition Probability','Entropy Rate', ...
               'Dwell Time','pMEM R2'};
cols = {'Metric','Analysis','Participant_ID','Band','Comparison','N_Obs','N_Participants', ...
    'Estimate','SE','CI_Lower','CI_Upper','Statistic','Statistic_Type', ...
    'DF1','DF2','P_Raw','P_FDR','FDR_Family','Model','Notes', ...
    'Standardized_Beta','Standardized_SE','Standardized_CI_Lower','Standardized_CI_Upper'};
rows = cell(0,numel(cols));

T.condition = categorical(T.condition,{'A','B'});
[G,~] = findgroups(T.subject);
subjMeanAcc = splitapply(@mean,T.accuracy,G);
T.accuracy_c = T.accuracy - subjMeanAcc(G);

for m = 1:5
    metric = metricNames{m};
    subT = T(T.measure==m,:);
    valid = isfinite(subT.value) & ~isundefined(subT.condition) & ...
            isfinite(subT.subject);
    if m ~= 5
        valid = valid & isfinite(subT.accuracy_c);
    end
    subT = subT(valid,:);
    if isempty(subT), error('No usable data for %s.',metric); end

    % Model specifications match the figure exactly.
    if m == 5
        model = 'value ~ condition + (1 | subject)';
    else
        model = 'value ~ condition + accuracy_c + (1 | subject)';
    end
    mdl = fitlme(subT,model);
    c = mdl.Coefficients;
    nObs = height(subT);
    nSub = numel(unique(subT.subject));
    outcomeSD = withinSD(subT.value,subT.subject);

    if m == 5
        modelLabel = 'pMEM R2 stimulation model';
    else
        modelLabel = 'Accuracy-adjusted stimulation LME';
    end

    ixCondition = find(strcmp(c.Name,'condition_B'),1);
    if isempty(ixCondition), error('Missing condition_B coefficient for %s.',metric); end
    rows(end+1,:) = coefRow(metric,modelLabel,'All participants', ...
        'ON - OFF','condition_B',subT,mdl,ixCondition,1/outcomeSD, ...
        'Condition difference in pooled within-patient outcome-SD units; unadjusted p-value.',model); %#ok<AGROW>

    % The accuracy covariate is included only in the four metric models.
    if m ~= 5
        ixAcc = find(strcmp(c.Name,'accuracy_c'),1);
        if isempty(ixAcc), error('Missing accuracy_c coefficient for %s.',metric); end
        accSD = withinSD(subT.accuracy_c,subT.subject);
        rows(end+1,:) = coefRow(metric,modelLabel,'All participants', ...
            'Accuracy covariate','accuracy_c',subT,mdl,ixAcc,accSD/outcomeSD, ...
            'Continuous predictor standardized with pooled within-patient SDs.',model); %#ok<AGROW>
    end

    % Descriptive comparisons reproduce figure z-scoring within patient
    % over the OFF and ON windows combined. Group SEM uses patients (n=2).
    groupZ = nan(nPatients,2);
    for s = 1:nPatients
        ptT = subT(subT.subject==s,:);
        if isempty(ptT), continue; end
        mu = mean(ptT.value);
        sigma = std(ptT.value);
        if sigma == 0, sigma = 1; end
        for cond = 1:2
            if cond == 1
                label = 'OFF'; ix = ptT.condition == 'A';
            else
                label = 'ON'; ix = ptT.condition == 'B';
            end
            vals = ptT.value(ix);
            if isempty(vals), continue; end
            z = (vals-mu)./sigma;
            groupZ(s,cond) = mean(z);
            rawSE = std(vals)/sqrt(numel(vals));
            rows(end+1,:) = makeRow(metric,'Patient ON/OFF descriptive',alllabels{s}, ...
                'beta',label,numel(vals),1,mean(vals),rawSE,NaN,NaN, ...
                NaN,'descriptive',NaN,NaN,NaN,NaN,'None','Raw condition summary', ...
                sprintf('Raw mean; patient z mean = %.8g; patient z SD = %.8g; window-level SEM is descriptive only.', ...
                mean(z),std(z))); %#ok<AGROW>
        end
    end
    for cond = 1:2
        if cond == 1, label = 'OFF'; else, label = 'ON'; end
        x = groupZ(:,cond);
        x = x(isfinite(x));
        if isempty(x), continue; end
        rows(end+1,:) = makeRow(metric,'Across-patient z summary','All participants', ...
            'beta',label,NaN,numel(x),mean(x),std(x)/sqrt(numel(x)),NaN,NaN, ...
            NaN,'descriptive',NaN,NaN,NaN,NaN,'None','Within-patient z means', ...
            'Mean and SEM of patient-specific ON/OFF z means (not window-weighted).'); %#ok<AGROW>
    end
end

results = cell2table(rows,'VariableNames',cols);
writetable(results,outFile,'Sheet','Figure4_Statistics','WriteMode','overwritesheet');
fprintf('Exported %d Figure 4 rows to %s\n',height(results),outFile);
end

function row = coefRow(metric,analysis,participant,comparison,coefName,t,mdl,ix,scale,note,model)
c = mdl.Coefficients;
if ~strcmp(c.Name{ix},coefName)
    error('Unexpected coefficient name for %s.',metric);
end
df = c.DF(ix);
crit = tinv(0.975,df);
est = c.Estimate(ix); se = c.SE(ix);
lo = est-crit*se; hi = est+crit*se;
row = makeRow(metric,analysis,participant,'beta',comparison,height(t), ...
    numel(unique(t.subject)),est,se,lo,hi,c.tStat(ix),'t',df,NaN,c.pValue(ix), ...
    NaN,'None',model,note,est*scale,se*scale,lo*scale,hi*scale);
end

function sd = withinSD(y,p)
num = 0; den = 0;
for id = unique(p(:))'
    v = y(p==id);
    v = v(isfinite(v));
    if numel(v) >= 2
        num = num + (numel(v)-1)*var(v);
        den = den + numel(v)-1;
    end
end
if den == 0 || num <= 0
    sd = NaN;
else
    sd = sqrt(num/den);
end
end

function row = makeRow(varargin)
row = varargin;
if numel(row) == 20
    row = [row,{NaN,NaN,NaN,NaN}];
elseif numel(row) ~= 24
    error('Unexpected reporting row width %d.',numel(row));
end
end

function [results,spectra] = exportFigure6Statistics(T_PI,T_PI_PCA,T_PCA, ...
    allData,pcaCutoff,alllabels,outFile)
% EXPORTFIGURE6STATISTICS  Detailed Figure 6 predictive-information report.
%
% Call after the three analysis tables and pcaCutoff exist in Fig6.m:
%   exportFigure6Statistics(T_PI,T_PI_PCA,T_PCA,allData,pcaCutoff, ...
%       alllabels,fullfile(outDir,'Figure6_Statistics.xlsx'));
%
% Sheets:
%   Figure6_Statistics  Three exact LMEs, coefficients, subject correlations,
%                       pooled visual correlations, cohort/screening details.
%   Figure6_PC_Spectra Per-subject least/most-depressed recordings, PC 1:20,
%                       and across-subject mean / SEM (as in the figure).
%
% LME fits use raw data. Fixed effects are standardized afterward using
% predictor and outcome pooled within-patient SDs from that model's dataset.
% The percentage attenuation is descriptive, NOT a mediation estimate.

if nargin < 7 || isempty(outFile)
    outFile = fullfile(pwd,'Figure6_Statistics.xlsx');
end
if exist(outFile,'file'), delete(outFile); end

cols = {'Metric','Analysis','Participant_ID','Band','Comparison','N_Obs','N_Participants', ...
    'Estimate','SE','CI_Lower','CI_Upper','Statistic','Statistic_Type', ...
    'DF1','DF2','P_Raw','P_FDR','FDR_Family','Model','Notes', ...
    'Standardized_Beta','Standardized_SE','Standardized_CI_Lower','Standardized_CI_Upper'};
rows = cell(0,numel(cols));
subjects = categories(allData.patient);
if numel(alllabels) ~= numel(subjects)
    error('Participant labels do not match the participants in allData.');
end

%% Fit the three models exactly as specified in Fig6.m
modelNames = {'Model 1: PI ~ severity', ...
              'Model 2: PI ~ severity + PCA', ...
              'Model 3: PCA ~ severity'};
modelForms = {'rawPI_max ~ cat + (1|patient)', ...
              'rawPI_max ~ cat + pcaVariance + (1|patient)', ...
              'pcaVariance ~ cat + (1|patient)'};
modelData = {T_PI,T_PI_PCA,T_PCA};
modelOutcomes = {'rawPI_max','rawPI_max','pcaVariance'};
modelMetrics = {'Raw predictive information','Raw predictive information','PCA variance explained'};
mdl = cell(3,1);
stdSeverity = nan(3,1);

for k = 1:3
    tbl = modelData{k};
    mdl{k} = fitlme(tbl,modelForms{k});
    c = mdl{k}.Coefficients;
    y = tbl.(modelOutcomes{k});
    ySD = withinSD(y,tbl.patient);
    for j = 1:height(c)
        coefficientName = c.Name{j};
        if strcmp(coefficientName,'cat')
            xSD = withinSD(tbl.cat,tbl.patient);
            factor = xSD/ySD;
        elseif strcmp(coefficientName,'pcaVariance')
            xSD = withinSD(tbl.pcaVariance,tbl.patient);
            factor = xSD/ySD;
        else
            factor = NaN;  % Do not standardize the intercept.
        end
        crit = tinv(0.975,c.DF(j));
        lo = c.Estimate(j)-crit*c.SE(j);
        hi = c.Estimate(j)+crit*c.SE(j);
        rows(end+1,:) = makeRow(modelMetrics{k},modelNames{k}, ...
            'All participants','Resting',coefficientName,height(tbl), ...
            numel(unique(tbl.patient)),c.Estimate(j),c.SE(j),lo,hi, ...
            c.tStat(j),'t',c.DF(j),NaN,c.pValue(j),NaN,'None', ...
            modelForms{k},sprintf('Raw-data LME; pcaVariance <= %.8g; standardized by model-specific pooled within-patient SDs.',pcaCutoff), ...
            c.Estimate(j)*factor,c.SE(j)*factor,lo*factor,hi*factor); %#ok<AGROW>
        if strcmp(coefficientName,'cat')
            stdSeverity(k) = c.Estimate(j)*factor;
        end
    end
end

%% Descriptive attenuation, matching the figure definition
if isfinite(stdSeverity(1)) && stdSeverity(1) ~= 0
    attenuation = 100*(stdSeverity(1)-stdSeverity(2))/stdSeverity(1);
else
    attenuation = NaN;
end
rows(end+1,:) = makeRow('Raw predictive information','Severity coefficient attenuation', ...
    'All participants','Resting','Model 1 versus Model 2',NaN,NaN, ...
    attenuation,NaN,NaN,NaN,NaN,'percent',NaN,NaN,NaN,NaN,'None', ...
    '100 * (betaStdModel1-betaStdModel2)/betaStdModel1', ...
    'Descriptive comparison of model-specific standardized coefficients; NOT a mediation estimate.');

%% Per-participant raw Pearson correlations and plotted pooled z correlations
corrDefinitions = {T_PI,'rawPI_max','Raw predictive information'; ...
                   T_PCA,'pcaVariance','PCA variance explained'};
for q = 1:size(corrDefinitions,1)
    tbl = corrDefinitions{q,1};
    feature = corrDefinitions{q,2};
    label = corrDefinitions{q,3};
    Xz = []; Yz = []; nIncluded = 0;
    for si = 1:numel(subjects)
        idx = tbl.patient == subjects{si};
        x = tbl.cat(idx); y = tbl.(feature)(idx);
        valid = isfinite(x) & isfinite(y);
        x = x(valid); y = y(valid);
        if numel(x)<2 || std(x)<=0 || std(y)<=0, continue; end
        if numel(x)>2
            [rv,pv] = corr(x,y,'Type','Pearson');
            rows(end+1,:) = makeRow(label,'Participant Pearson correlation', ...
                alllabels{si},'Resting','Severity vs metric',numel(x),1, ...
                rv,NaN,NaN,NaN,NaN,'Pearson r',NaN,NaN,pv,NaN,'None', ...
                'Pearson correlation','Raw recordings; no multiple-testing correction in figure.'); %#ok<AGROW>
        end
        nIncluded = nIncluded + 1;
        Xz = [Xz;(x-mean(x))/std(x)]; %#ok<AGROW>
        Yz = [Yz;(y-mean(y))/std(y)]; %#ok<AGROW>
    end
    if ~isempty(Xz)
        [rv,pv] = corr(Xz,Yz,'Type','Pearson');
        rows(end+1,:) = makeRow(label,'Pooled visual correlation', ...
            'All participants','Resting','Within-patient z scores',numel(Xz), ...
            nIncluded,rv,NaN,NaN,NaN,NaN,'Pearson r',NaN,NaN, ...
            pv,NaN,'None','Pearson correlation of within-patient z scores', ...
            'Descriptive pooled correlation shown in figure; not a hierarchical LME p-value.'); %#ok<AGROW>
    end
end

%% Sample counts, for auditability of the three actual model datasets
for k = 1:3
    tbl = modelData{k};
    rows(end+1,:) = makeRow(modelMetrics{k},'Model sample count','All participants', ...
        'Resting',modelNames{k},height(tbl),numel(unique(tbl.patient)), ...
        height(tbl),NaN,NaN,NaN,NaN,'N',NaN,NaN,NaN,NaN,'None', ...
        modelForms{k},sprintf('Complete-case observations after global 90th-percentile PCA cutoff %.8g.',pcaCutoff)); %#ok<AGROW>
    for si = 1:numel(subjects)
        n = nnz(tbl.patient == subjects{si});
        rows(end+1,:) = makeRow(modelMetrics{k},'Participant sample count', ...
            alllabels{si},'Resting',modelNames{k},n,1,n,NaN,NaN,NaN, ...
            NaN,'N',NaN,NaN,NaN,NaN,'None',modelForms{k}, ...
            'Participant observations in exact model dataset.'); %#ok<AGROW>
    end
end

results = cell2table(rows,'VariableNames',cols);
writetable(results,outFile,'Sheet','Figure6_Statistics','WriteMode','overwritesheet');

%% PCA variance spectra: patient extrema after the same cutoff
% Reflects the tile that averages PC 1:20 over the six participants.
nPC = 20;
spectraRows = cell(0,7);
varLeast = nan(numel(subjects),nPC);
varMost = nan(numel(subjects),nPC);
for si = 1:numel(subjects)
    t = allData(allData.patient==subjects{si},:);
    t = t(isfinite(t.cat),:);
    if isempty(t), continue; end
    [catLeast,iLeast] = min(t.cat);
    [catMost,iMost] = max(t.cat);
    for which = 1:2
        if which==1
            v = t.VarExp{iLeast}; severity = catLeast;
            condition = 'Least depressed';
        else
            v = t.VarExp{iMost}; severity = catMost;
            condition = 'Most depressed';
        end
        v = v(:);
        for pc = 1:min(nPC,numel(v))
            spectraRows(end+1,:) = {alllabels{si},condition,pc,severity,v(pc),NaN,1}; %#ok<AGROW>
        end
        if which == 1
            varLeast(si,1:min(nPC,numel(v))) = v(1:min(nPC,numel(v)));
        else
            varMost(si,1:min(nPC,numel(v))) = v(1:min(nPC,numel(v)));
        end
    end
end
for which = 1:2
    if which == 1
        vmat = varLeast; condition = 'Least depressed';
    else
        vmat = varMost; condition = 'Most depressed';
    end
    for pc = 1:nPC
        v = vmat(:,pc);
        v = v(isfinite(v));
        if isempty(v), continue; end
        spectraRows(end+1,:) = {'All participants',condition,pc,NaN, ...
            mean(v),std(v)/sqrt(numel(v)),numel(v)}; %#ok<AGROW>
    end
end
spectra = cell2table(spectraRows,'VariableNames', ...
    {'Participant_ID','Severity_Group','PC','CAT_DI','Variance_Explained_Percent','SEM','N_Participants'});
writetable(spectra,outFile,'Sheet','Figure6_PC_Spectra','WriteMode','overwritesheet');
fprintf('Exported %d Figure 6 statistical rows and %d PC spectrum rows to %s\n', ...
    height(results),height(spectra),outFile);
end

function sd = withinSD(y,p)
num = 0; den = 0;
for id = categories(p)'
    v = y(p==id{1});
    v = v(isfinite(v));
    if numel(v)>=2
        num = num+(numel(v)-1)*var(v);
        den = den+numel(v)-1;
    end
end
if den==0
    sd = NaN;
else
    sd = sqrt(num/den);
end
end

function row = makeRow(varargin)
row = varargin;
if numel(row)==20
    row = [row,{NaN,NaN,NaN,NaN}];
elseif numel(row)~=24
    error('Unexpected reporting row width: %d',numel(row));
end
end

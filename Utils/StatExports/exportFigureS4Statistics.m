function results = exportFigureS4Statistics(T, alllabels, outFile)
% EXPORTFIGURES4STATISTICS  Detailed statistics for Figure S4 (AmpModControl).
%
% Re-fits exactly the single-participant regressions in the revised FigS4:
%   Measures 1-4: value ~ condition + accuracy_c
%   Measure 5:   value ~ condition               (pMEM R^2)
%
% INPUTS
%   T         Long-form table built by FigS4, with columns:
%             subject, measure, condition, trial, value, accuracy, accuracy_c.
%             Condition A = OFF; condition B = ON.
%   alllabels Participant display labels from FigS4 (one patient).
%   outFile   Destination .xlsx file.
%
% OUTPUT
%   results   One rectangular table on the FigureS4_Statistics worksheet.
%             Includes primary ON-minus-OFF effects, remaining coefficients,
%             model details, sample counts, means, SEs, 95% CIs and p-values.
%
% No FDR adjustment is used in FigS4, so P_FDR is left NaN; no new
% multiplicity adjustment is introduced by this exporter.
%
% NOTE: The p-values treat observation windows as independent, matching the
% figure's ordinary linear regressions; they do not adjust for serial
% correlation across adjacent windows.

if nargin < 3 || isempty(outFile)
    outFile = fullfile(pwd, 'FigureS4_Statistics.xlsx');
end

required = {'subject','measure','condition','value','accuracy','accuracy_c'};
assert(istable(T), 'T must be a table.');
assert(all(ismember(required, T.Properties.VariableNames)), ...
    'T is missing one or more required columns: %s', strjoin(required, ', '));
assert(numel(alllabels) == 1, ...
    'Figure S4 is a single-participant analysis: expected one display label.');
assert(numel(unique(T.subject)) == 1, ...
    'Figure S4 must have observations from exactly one participant.');

% Preserve FigS4's reference category A (OFF), so the condition_B
% coefficient is the adjusted ON-minus-OFF contrast.
T.condition = categorical(T.condition, {'A','B'});
assert(all(~ismissing(T.condition)), ...
    'Condition must be coded as A (OFF) or B (ON).');

metricNames = {'Temperature','Transition Probability', ...
               'Entropy Rate','Dwell Time','pMEM R2'};
cols = {'Metric','Analysis','Participant_ID','Band','Comparison', ...
    'N_Obs','N_OFF','N_ON','Estimate','SE','CI_Lower','CI_Upper', ...
    'Statistic','Statistic_Type','DF1','DF2','P_Raw','P_FDR', ...
    'FDR_Family','Model','Notes','R2','Adjusted_R2','RMSE', ...
    'OFF_Mean','ON_Mean'};
rows = cell(0, numel(cols));

for m = 1:5
    subT = T(T.measure == m, :);

    if m == 5
        formula = 'value ~ condition';
        modelText = 'value ~ condition';
        valid = isfinite(subT.value) & ~ismissing(subT.condition);
    else
        formula = 'value ~ condition + accuracy_c';
        modelText = 'value ~ condition + accuracy_c';
        valid = isfinite(subT.value) & isfinite(subT.accuracy_c) & ...
            ~ismissing(subT.condition);
    end
    % fitlm excludes incomplete rows. Matching this subset makes counts and
    % summary means refer to the same observations as the fitted model.
    subValid = subT(valid, :);
    nOff = nnz(subValid.condition == 'A');
    nOn  = nnz(subValid.condition == 'B');
    assert(nOff > 0 && nOn > 0, ...
        'Measure %d must contain valid observations in both conditions.', m);

    lm = fitlm(subT, formula);
    c = lm.Coefficients;
    names = lm.CoefficientNames;
    if isstring(names), names = cellstr(names); end
    assert(height(c) == numel(names), ...
        'Coefficient name count does not match table for measure %d.', m);
    assert(lm.NumObservations == height(subValid), ...
        'Valid-row count does not match fitted sample size for measure %d.', m);

    idxCond = find(strcmp(names, 'condition_B'), 1);
    assert(~isempty(idxCond), ...
        'condition_B coefficient missing for measure %d. Check condition coding.', m);

    offMean = mean(subValid.value(subValid.condition == 'A'));
    onMean  = mean(subValid.value(subValid.condition == 'B'));
    tc = tinv(0.975, lm.DFE);
    r2 = lm.Rsquared.Ordinary;
    r2adj = lm.Rsquared.Adjusted;
    rmse = lm.RMSE;

    for j = [idxCond, setdiff(1:numel(names),idxCond,'stable')]
        est = c.Estimate(j);
        se = c.SE(j);
        if j == idxCond
            analysis = 'Primary ON vs OFF';
            comparison = 'ON - OFF (B - A)';
            if m == 5
                note = 'Unadjusted condition effect; R2 is the dependent variable.';
            else
                note = 'Condition effect adjusted for within-participant centered pMEM R2.';
            end
        else
            analysis = 'Other model coefficient';
            comparison = names{j};
            note = 'Additional coefficient from exactly the same fitted model.';
        end

        rows(end+1,:) = {metricNames{m},analysis,alllabels{1},'beta', ...
            comparison,height(subValid),nOff,nOn,est,se, ...
            est-tc*se,est+tc*se,c.tStat(j),'t',1,lm.DFE, ...
            c.pValue(j),NaN,'None',modelText,note, ...
            r2,r2adj,rmse,offMean,onMean}; %#ok<AGROW>
    end
end

results = cell2table(rows, 'VariableNames', cols);
[outDir,~,~] = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
if exist(outFile,'file'), delete(outFile); end
writetable(results, outFile, 'Sheet', 'FigureS4_Statistics');
fprintf('Figure S4: exported %d rows to %s\n', height(results), outFile);
end

%% FIGURE S4 | Single-participant amplitude-modulation control
% AmpModControl remains the input dataset. Models use ordinary
% regression instead of a random intercept with one participant.

% SINGLE-PATIENT CONTROL ANALYSIS

allpts = {''};
alllabels = {''};
alltasks = {'AmpModControl'};
cmap = lines(5);

allbands = {'beta'};
thresholds = round([50]);

% FIGURE: timeline on left, 2x2 metric grid on right
% 6 x 8 grid:
%   left 4 columns = timeline
%   right 4 columns = 2x2 metric plots

figure( ...
    'Renderer','painters', ...
    'Position',[10 10 1000 300], ...
    'Color','w');

outerTL = tiledlayout(6,8, ...
    'TileSpacing','compact', ...
    'Padding','compact');

stimonvals = {};
stimoffvals = {};

% LOAD SINGLE PATIENT

pt = 1;

valson = [];
valsoff = [];

for thresh = 1:length(thresholds)

    for tk = 1:length(alltasks)

        for b = 1:length(allbands)

            paramfile = betterLoadSave([ ...
                USER_DIR,'/Analysis/', ...
                alltasks{tk}, '/', ...
                allpts{pt}, ...
                '/Params/pMEM_features_', ...
                num2str(thresholds(thresh)), '/', ...
                allbands{b}, '_params.mat']);

            load(paramfile);

            % Sort windows chronologically

            segnum = cellfun(@str2num, ...
                extractBetween({params.filename},'Window_','O'));

            [~,sortorder] = sort(segnum);

            paramsordered = params(sortorder);

            stimval = double( ...
                contains({paramsordered.filename},'ON'));

            % Extract metrics

            paramsParameters = [paramsordered.parameters];

            temps = [paramsParameters.effectiveTemperature];

            mets = [paramsordered.metrics];
            mets = [mets.sharedMetrics];

            R2pmem = [mets.R2pmem];

            Rall = zeros(1,length(paramsordered));
            meanDwellAll = zeros(1,length(paramsordered));
            HrateAll = zeros(1,length(paramsordered));

            ticklocs = zeros(1,length(paramsordered));
            winlocs = zeros(length(paramsordered),2);

            maxt = [];

            % Extract dynamics and timing

            for pf = 1:length(paramsordered)

                Dtrans = ...
                    paramsordered(pf).sharedFeatures.Dynamics.Dtrans;

                [ ...
                    Rall(pf), ...
                    meanDwellAll(pf), ...
                    HrateAll(pf), ...
                    ~] = summarizeTransitions(Dtrans);

                load( ...
                    betterLoadSave(paramsordered(pf).filename), ...
                    'timevec');

                ticklocs(pf) = mean(timevec);

                winlocs(pf,1) = timevec(1);
                winlocs(pf,2) = timevec(end);

                if pf == length(paramsordered)
                    maxt = timevec(end);
                end

            end

            ticklocs = ticklocs/60;
            maxt = maxt/60;
            winlocs = winlocs/60;

            % Normalize metrics for timeline visualization

            temp_n = (temps - min(temps)) ./ ...
                (max(temps)-min(temps));

            tp_n = (Rall - min(Rall)) ./ ...
                (max(Rall)-min(Rall));

            dt_n = (meanDwellAll - min(meanDwellAll)) ./ ...
                (max(meanDwellAll)-min(meanDwellAll));

            er_n = (HrateAll - min(HrateAll)) ./ ...
                (max(HrateAll)-min(HrateAll));

            temp_plot = temp_n;
            tp_plot   = tp_n + 1.5;
            dt_plot   = dt_n + 3.0;
            er_plot   = er_n + 4.5;

            % TIMELINE
            % Single patient occupies 2 x 4 tiles on left.

            axTop = nexttile( ...
                outerTL, ...
                1, ...
                [6 4]);

            hold(axTop,'on');

            % Highlight stimulation ON windows

            d = diff([0 stimval 0]);

            startIdx = find(d == 1);
            endIdx   = find(d == -1)-1;

            yl = [-0.2 5.7];

            for k = 1:numel(startIdx)

                x1 = winlocs(startIdx(k),1);
                x2 = winlocs(endIdx(k),2);

                patch( ...
                    axTop, ...
                    [x1 x2 x2 x1], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    [1 0.85 0.85], ...
                    'EdgeColor','none', ...
                    'FaceAlpha',0.4);

            end

            % Plot normalized metrics

            yyaxis(axTop,'left');

            plot( ...
                axTop, ...
                ticklocs,temp_plot, ...
                '-', ...
                'Color',cmap(1,:), ...
                'LineWidth',2);

            plot( ...
                axTop, ...
                ticklocs,tp_plot, ...
                '-', ...
                'Color',cmap(2,:), ...
                'LineWidth',2);

            plot( ...
                axTop, ...
                ticklocs,dt_plot, ...
                '-', ...
                'Color',cmap(5,:), ...
                'LineWidth',2);

            plot( ...
                axTop, ...
                ticklocs,er_plot, ...
                '-', ...
                'Color',cmap(3,:), ...
                'LineWidth',2);

            ylim(axTop,yl);
            xlim(axTop,[0 maxt]);

            title(axTop,alllabels{pt});
            xlabel(axTop,'Task Time (Minute)');

            box(axTop,'off');

            % White separators

            line(axTop,[0 0],[2.505 2.995], ...
                'Color','w','LineWidth',1.5);

            line(axTop,[0 0],[1.005 1.495], ...
                'Color','w','LineWidth',1.5);

            line(axTop,[0 0],[-0.1995 -0.005], ...
                'Color','w','LineWidth',1.5);

            line(axTop,[0 0],[4.005 4.495], ...
                'Color','w','LineWidth',1.5);

            line(axTop,[0 0],[5.505 5.7], ...
                'Color','w','LineWidth',1.5);

            % Left-axis labels

            yticks(axTop,[0 1 1.5 2.5 3 4 4.5 5.5]);

            temp_min = min(temps);
            temp_max = max(temps);

            tp_min = min(Rall);
            tp_max = max(Rall);

            dt_min = min(meanDwellAll);
            dt_max = max(meanDwellAll);

            ent_min = min(HrateAll);
            ent_max = max(HrateAll);

            yticklabels(axTop,{ ...
                num2str(temp_min,'%.3g')
                num2str(temp_max,'%.3g')
                num2str(tp_min,'%.3g')
                num2str(tp_max,'%.3g')
                num2str(dt_min,'%.3g')
                num2str(dt_max,'%.3g')
                num2str(ent_min,'%.3g')
                num2str(ent_max,'%.3g')});

            ylim(axTop,[-0.2 5.7]);

            % pMEM R2

            yyaxis(axTop,'right');

            plot( ...
                axTop, ...
                ticklocs,R2pmem, ...
                'k', ...
                'LineWidth',2);

            ylabel(axTop,'pMEM R^2');
            ylim(axTop,[0 1]);

            axTop.YAxis(1).Color = 'k';
            axTop.YAxis(2).Color = 'k';

            % Store values

            valson(1,:) = temps(stimval==1);
            valsoff(1,:) = temps(stimval==0);

            valson(2,:) = Rall(stimval==1);
            valsoff(2,:) = Rall(stimval==0);

            valson(3,:) = HrateAll(stimval==1);
            valsoff(3,:) = HrateAll(stimval==0);

            valson(4,:) = meanDwellAll(stimval==1);
            valsoff(4,:) = meanDwellAll(stimval==0);

            valson(5,:) = R2pmem(stimval==1);
            valsoff(5,:) = R2pmem(stimval==0);

        end
    end
end

stimonvals{pt} = valson;
stimoffvals{pt} = valsoff;

% BUILD LONG-FORM TABLE

T = table();

A = stimoffvals{1};
B = stimonvals{1};

Ni = size(A,2);
Nj = size(B,2);

for m = 1:5

    TA = table( ...
        repmat(1,Ni,1), ...
        repmat(m,Ni,1), ...
        categorical(repmat("A",Ni,1)), ...
        (1:Ni)', ...
        A(m,:)', ...
        A(5,:)');

    TB = table( ...
        repmat(1,Nj,1), ...
        repmat(m,Nj,1), ...
        categorical(repmat("B",Nj,1)), ...
        (1:Nj)', ...
        B(m,:)', ...
        B(5,:)');

    T = [T; TA; TB];

end

T.Properties.VariableNames = ...
    {'subject','measure','condition','trial','value','accuracy'};

disp('Long-form table T successfully built.');

% R2 MODEL (NO PATIENT RANDOM EFFECT)

T.condition = categorical( ...
    T.condition, ...
    {'A','B'});

[G,~] = findgroups(T.subject);

subjMeanAcc = splitapply( ...
    @mean, ...
    T.accuracy, ...
    G);

% Center R2 across the single participant's windows; the same
% centered covariate is used in each of the four metric models.
T.accuracy_c = ...
    T.accuracy - subjMeanAcc(G);

subT_acc = T(T.measure==5,:);

lm_R2pmem = fitlm( ...
    subT_acc, ...
    'value ~ condition');

disp('--------------------------------------');
disp('R2 MODEL (SINGLE PATIENT)');
disp(lm_R2pmem);

% ACCURACY-ADJUSTED MODELS (SINGLE-PATIENT REGRESSION)


% Fit accuracy-adjusted linear models for the four landscape metrics.
% Single-patient analysis: no subject random effect.

metricModels = cell(1,4);
pvals_condition = nan(1,4);

for m = 1:4

    subT = T(T.measure == m, :);

    % Model the effect of stimulation while adjusting for accuracy.
    lm = fitlm(subT, 'value ~ condition + accuracy_c');

    metricModels{m} = lm;

    % Extract the ON versus OFF stimulation coefficient.
    coefTable = lm.Coefficients;
    condIdx = strcmp(lm.CoefficientNames, 'condition_B');

    assert(any(condIdx), ...
        'Could not find condition_B coefficient for measure %d.', m);

    pvals_condition(m) = coefTable.pValue(condIdx);

    % Display model results.
    disp('--------------------------------------');
    fprintf('MEASURE %d MODEL\n', m);
    disp(lm);

end


% METRIC PLOTS
% 2 x 2 grid occupying the right half.

measurelabels = { ...
    'Temperature', ...
    'Transition Probability', ...
    'Entropy Rate', ...
    'Dwell Time'};

% 6 x 8 layout:
% row 1-3: metric 1 | metric 2
% row 4-6: metric 3 | metric 4
% Each metric occupies 3 x 2 tiles.

metricStarts = [5 7 29 31];

for m = 1:4

    axBot = nexttile( ...
        outerTL, ...
        metricStarts(m), ...
        [3 2]);

    hold(axBot,'on');

    subT = T(T.measure==m,:);

    % Z-score within subject

    sT = subT(subT.subject==1,:);

    mu = mean(sT.value);
    sigma = std(sT.value);

    if sigma == 0
        sigma = 1;
    end

    Avals = ...
        (sT.value(sT.condition=='A') - mu) ./ sigma;

    Bvals = ...
        (sT.value(sT.condition=='B') - mu) ./ sigma;

    % Individual observations

    scatter( ...
        axBot, ...
        1 + 0.05*randn(size(Avals)), ...
        Avals, ...
        30, ...
        'MarkerFaceColor',cmap(3,:), ...
        'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.65);

    scatter( ...
        axBot, ...
        2 + 0.05*randn(size(Bvals)), ...
        Bvals, ...
        30, ...
        'MarkerFaceColor',cmap(3,:), ...
        'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.65);

    % Subject mean

    meanA = mean(Avals);
    meanB = mean(Bvals);

    plot( ...
        axBot, ...
        [1 2], ...
        [meanA meanB], ...
        '--', ...
        'Color',cmap(3,:), ...
        'LineWidth',1.5);

    % Mean line

    plot( ...
        axBot, ...
        [1 2], ...
        [meanA meanB], ...
        'k-', ...
        'LineWidth',2);

    errorbar( ...
        axBot, ...
        1, ...
        meanA, ...
        0, ...
        'k', ...
        'LineWidth',1.5);

    errorbar( ...
        axBot, ...
        2, ...
        meanB, ...
        0, ...
        'k', ...
        'LineWidth',1.5);

    % p-value

    yl = ylim(axBot);

    txt = sprintf( ...
        'p = %.3f', ...
        pvals_condition(m));

    text( ...
        axBot, ...
        1.05, ...
        yl(2)-0.05*(yl(2)-yl(1)), ...
        txt, ...
        'FontSize',10);

    % Formatting

    set(axBot, ...
        'XTick',[1 2], ...
        'XTickLabel',{'OFF','ON'});

    xlabel(axBot,'Stimulation');
    ylabel(axBot,[measurelabels{m} ' (z)']);

    xlim(axBot,[0.8 2.2]);

    box(axBot,'off');

end

% SAVE FIGURE S4
statsDir = [USER_DIR,'/StatTablesNew'];
S4stats = exportFigureS4Statistics( ...
    T, alllabels, ...
    fullfile(statsDir, 'FigureS4_Statistics.xlsx'));
exportgraphics(gcf, betterLoadSave( [USER_DIR,'/PaperFigures/FigS4_Stimulation.svg']));

close all;

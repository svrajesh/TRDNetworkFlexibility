
allpts = {};
alllabels = {};
alltasks = {'AmpMod'};
cmap = lines(5);

USER_DIR = '';

outDir = fullfile(rootDir,'PaperFigures');
allbands = {'beta'};
thresholds = round([50]);

%% Figure 4: Stimulation effects in two patients
% Two patients, beta-band pMEM metrics. Left: individual recording timelines.
% Right: four OFF-versus-ON comparisons, with individual and group summaries.
% Statistical models use raw metrics; displayed summaries are patient z-scores.

figure( ...
    'Renderer','painters', ...
    'Position',[10 10 1000 550], ...
    'Color','w');

% 4 x 8 grid: left four columns for two stacked timelines (2 rows each),
% right four columns for a 2 x 2 array of metric comparisons.
outerTL = tiledlayout(4,8, ...
    'TileSpacing','compact', ...
    'Padding','compact');

stimonvals = {};
stimoffvals = {};

%% Load stimulation windows and collect patient-level values
for pt = 1:length(allpts)

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

                % Order recording windows by their numeric filename index.
                segnum = cellfun(@str2num, ...
                    extractBetween({params.filename},'Window_','O'));

                [~,sortorder] = sort(segnum);
                paramsordered = params(sortorder);

                % Stimulation condition: ON = 1, OFF = 0.
                stimval = double( ...
                    contains({paramsordered.filename},'ON'));

                % Temperature and pMEM fit accuracy (R2) for each window.
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

                % Transition dynamics and original timing of each window.
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

                % Convert the timing information to minutes.
                ticklocs = ticklocs/60;
                maxt = maxt/60;
                winlocs = winlocs/60;

                % Normalize each metric to [0, 1] within this patient, then
                % vertically offset the curves for visualization only.
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

                % Each patient's timeline spans two rows by four columns.
                topRow = 1 + (pt-1)*2;

                axTop = nexttile( ...
                    outerTL, ...
                    (topRow-1)*8 + 1, ...
                    [2 4]);

                hold(axTop,'on');

                % Shade contiguous ON periods using actual window boundaries.
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

                % Four normalized and vertically separated metric traces.
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

                title( ...
                    axTop, ...
                    alllabels{pt});

                xlabel( ...
                    axTop, ...
                    'Task Time (Minute)');

                box(axTop,'off');

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

                % Show each metric's original minimum and maximum as tick labels.
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

                % Overlay pMEM R2 on a separate, unnormalized axis.
                yyaxis(axTop,'right');

                plot( ...
                    axTop, ...
                    ticklocs,R2pmem, ...
                    'k', ...
                    'LineWidth',2);

                ylabel( ...
                    axTop, ...
                    'pMEM R^2');

                ylim(axTop,[0 1]);

                axTop.YAxis(1).Color = 'k';
                axTop.YAxis(2).Color = 'k';

                % Store raw metric values by stimulation condition.
                % Rows: temperature, escape probability, entropy rate,
                % mean dwell time, and pMEM R2.
                valson(1,:) = temps(stimval==1);
                valsoff(1,:) = temps(stimval==0);

                valson(2,:) = Rall(stimval==1);
                valsoff(2,:) = Rall(stimval==0);

                valson(4,:) = meanDwellAll(stimval==1);
                valsoff(4,:) = meanDwellAll(stimval==0);

                valson(3,:) = HrateAll(stimval==1);
                valsoff(3,:) = HrateAll(stimval==0);

                valson(5,:) = R2pmem(stimval==1);
                valsoff(5,:) = R2pmem(stimval==0);

            end
        end
    end

    stimonvals{pt} = valson;
    stimoffvals{pt} = valsoff;

end

%% Assemble long-form table for mixed-effects models
% A = OFF, B = ON; one row per recording window and outcome metric.
% The accuracy column repeats the corresponding window's pMEM R2.
T = table();

for s = 1:length(allpts)

    A = stimoffvals{s};
    B = stimonvals{s};

    Ni = size(A,2);
    Nj = size(B,2);

    for m = 1:5

        TA = table( ...
            repmat(s,Ni,1), ...
            repmat(m,Ni,1), ...
            categorical(repmat("A",Ni,1)), ...
            (1:Ni)', ...
            A(m,:)', ...
            A(5,:)');

        TB = table( ...
            repmat(s,Nj,1), ...
            repmat(m,Nj,1), ...
            categorical(repmat("B",Nj,1)), ...
            (1:Nj)', ...
            B(m,:)', ...
            B(5,:)');

        T = [T; TA; TB];

    end
end

T.Properties.VariableNames = ...
    {'subject','measure','condition','trial','value','accuracy'};

disp('Long-form table T successfully built.');

%% Center pMEM fit accuracy within patient
% The centered R2 covariate is used only in the four metric models below.
T.condition = categorical( ...
    T.condition, ...
    {'A','B'});

[G,~] = findgroups(T.subject);

subjMeanAcc = splitapply( ...
    @mean, ...
    T.accuracy, ...
    G);

T.accuracy_c = ...
    T.accuracy - subjMeanAcc(G);

%% Model stimulation's association with pMEM R2 itself
% R2 is the outcome here; no accuracy covariate is included.
subT_acc = T(T.measure==5,:);

lme_R2pmem = fitlme( ...
    subT_acc, ...
    'value ~ condition + (1 | subject)');

disp('--------------------------------------');
disp('ACCURACY MODEL');
disp(lme_R2pmem);

%% Fit accuracy-adjusted stimulation models for the four metrics
% value ~ condition + accuracy_c + (1 | subject).
% The condition_B coefficient tests ON versus OFF at fixed pMEM R2.
% Keep the corresponding fit in lme_accuracy for unambiguous naming.
LMEmodels = cell(1,4);
pvals_condition = zeros(1,4);

for m = 1:4

    subT = T(T.measure==m,:);

    lme_accuracy = fitlme( ...
        subT, ...
        'value ~ condition + accuracy_c + (1 | subject)');

    LMEmodels{m} = lme_accuracy;

    % Store the ON-versus-OFF p-value from the adjusted model.
    coefTable = lme_accuracy.Coefficients;

    condIdx = strcmp( ...
        coefTable.Name, ...
        'condition_B');

    pvals_condition(m) = ...
        coefTable.pValue(condIdx);

    disp('--------------------------------------');
    fprintf('MEASURE %d MODEL\n',m);
    disp(lme_accuracy);

end

%% Plot OFF versus ON metric values for the two patients
% Four panels in the right half of the figure (2 x 2 layout).
measurelabels = { ...
    'Temperature', ...
    'Transition Probability', ...
    'Entropy Rate', ...
    'Dwell Time'};

subjColors = cmap(4:5,:);

% Upper-left, upper-right, lower-left, and lower-right panel starts.
metricStarts = [5 7 21 23];

for m = 1:4

    axBot = nexttile( ...
        outerTL, ...
        metricStarts(m), ...
        [2 2]);

    hold(axBot,'on');

    subT = T(T.measure==m,:);

    meanA = zeros(2,1);
    meanB = zeros(2,1);

    for s = 1:length(allpts)

        sT = subT(subT.subject==s,:);

        % Z-score across this patient's ON and OFF windows combined.
        mu = mean(sT.value);
        sigma = std(sT.value);

        if sigma == 0
            sigma = 1;
        end

        Avals = ...
            (sT.value(sT.condition=='A') - mu) ./ sigma;

        Bvals = ...
            (sT.value(sT.condition=='B') - mu) ./ sigma;

        % Individual window values, horizontally jittered for visibility.
        scatter( ...
            axBot, ...
            1 + 0.05*randn(size(Avals)), ...
            Avals, ...
            30, ...
            'MarkerFaceColor',subjColors(s,:), ...
            'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',0.65);

        scatter( ...
            axBot, ...
            2 + 0.05*randn(size(Bvals)), ...
            Bvals, ...
            30, ...
            'MarkerFaceColor',subjColors(s,:), ...
            'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',0.65);

        % Connect each patient's mean OFF and ON values.
        meanA(s) = mean(Avals);
        meanB(s) = mean(Bvals);

        plot( ...
            axBot, ...
            [1 2], ...
            [meanA(s) meanB(s)], ...
            '--', ...
            'Color',subjColors(s,:), ...
            'LineWidth',1.5);

    end

    % Group mean and SEM across the two patients (not windows).
    plot( ...
        axBot, ...
        [1 2], ...
        [mean(meanA) mean(meanB)], ...
        'k-', ...
        'LineWidth',2);

    errorbar( ...
        axBot, ...
        1, ...
        mean(meanA), ...
        std(meanA)/sqrt(2), ...
        'k', ...
        'LineWidth',1.5);

    errorbar( ...
        axBot, ...
        2, ...
        mean(meanB), ...
        std(meanB)/sqrt(2), ...
        'k', ...
        'LineWidth',1.5);

    yl = ylim(axBot);

    % Annotate the accuracy-adjusted stimulation p-value.
    txt = sprintf( ...
        'p = %.3e', ...
        pvals_condition(m));

    text( ...
        axBot, ...
        1.05, ...
        yl(2)-0.05*(yl(2)-yl(1)), ...
        txt, ...
        'FontSize',10);

    set(axBot, ...
        'XTick',[1 2], ...
        'XTickLabel',{'OFF','ON'});

    xlabel(axBot,'Stimulation');
    ylabel(axBot,[measurelabels{m} ' (z)']);

    xlim(axBot,[0.8 2.2]);

    box(axBot,'off');

end
exportFigure4Statistics(T,alllabels, ...
    fullfile(USER_DIR,'/StatTables', ...
             'Figure4_Statistics.xlsx'));

%% Export the two-patient stimulation figure
saveas(gcf,betterLoadSave([outDir,'Fig4_Stimulation.svg']));

close all;

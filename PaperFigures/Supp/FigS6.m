%% FIGURE S6 | Within-patient correlations of pMEM structure
% Original data selection, calculations, and inference retained.


%% Gather Info
allpts = {};
alllabels = {};
alltasks = {'Resting'};
allbands = {'beta'};
thresholds = [50];

for tk = 1:length(alltasks)
    task = alltasks{tk};
    for b = 1:length(allbands)
        band = allbands{b};
        for thresh = 1:length(thresholds)
            threshold = thresholds(thresh);
            figure('Renderer', 'painters', 'Position', [100, 100, 1000, 500], 'Color', 'w');
            tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            for pt = 1:length(allpts)
                pat = allpts{pt};
                paramsfile = betterLoadSave([USER_DIR,'/Analysis/',task,'/',pat,'/Params/pMEM_features_',num2str(threshold),'/',band,'_params.mat']);
                load(paramsfile);
                nodeNum = length(params(1).parameters.h);
                allhvals = zeros(length(params),nodeNum);
                allJvals = zeros(length(params),nodeNum,nodeNum);

                %% Gather the data
                for rec = 1:length(params)

                    allhvals(rec, :) = params(rec).parameters.h; % Store h values
                    allJvals(rec, :, :) = params(rec).parameters.J; % Store J values
                end
                pairwiseStructureCorrs = structureCorrelations(allhvals,allJvals);
                mask = triu(true(size(pairwiseStructureCorrs)),1);

                nexttile;
                % Use the mask to extract upper triangle elements
                upperTri = pairwiseStructureCorrs(mask);

                histogram(upperTri(:), 'Normalization', 'probability'); % Plot histogram

                % Calculate and plot the mean correlation coefficient
                meanCorr = mean(upperTri);
                hold on; % Hold the current plot
                yLimits = ylim; % Get current y-axis limits
                plot([meanCorr, meanCorr], yLimits, 'k:', 'LineWidth', lnWidth); % Add black dotted line
                hold off; % Release the hold on the current plot

                % Add text to indicate the mean correlation coefficient
                text(meanCorr - 0.08, yLimits(2) * 0.9, sprintf('$\\mu$: %.2f', meanCorr), 'Color', 'k', 'FontSize', fontSize, 'HorizontalAlignment', 'right','Interpreter','latex');
                xlim([-1,1]);
                xlabel('Correlation Coefficient');
                ylabel('Probability');
                title(alllabels{pt});
                box off

            end
            sgtitle('Histogram of Pairwise Structure Correlations');

           exportgraphics(gcf,betterLoadSave([USER_DIR,'/PaperFigures/FigS6_StructureCorrsRest.svg']), 'Resolution', 600);
        end
    end
end

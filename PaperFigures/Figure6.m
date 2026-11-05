%% =========================================================
% FIGURE 6: Predictive Information + PCA Variance Explained
%
% Primary metric:
%   rawPI_max
%
% Questions:
%
%   1. Does greater severity change PImax?
%
%   2. Does greater severity change PCA variance?
%
%   3. Does PCA variance account for the severity-related
%      change in PImax, or does PImax retain an independent
%      association with severity?
%
% Models:
%
%   Model 1:
%       rawPI_max ~ cat + (1|patient)
%
%   Model 2:
%       rawPI_max ~ cat + pcaVariance + (1|patient)
%
%   Model 3:
%       pcaVariance ~ cat + (1|patient)
%
% LME models are fit on RAW data.
% Standardized effects are calculated afterward using
% pooled within-patient SDs from the exact dataset used
% for each model.
%% =========================================================


%% =========================================================
% Load all predictive information results
%% =========================================================

allpts = {};

alllabels = {};



USER_DIR = '';
resultDir = betterLoadSave([USER_DIR,'Analysis/Predictability/Resting']);
outDir = fullfile(rootDir,'PaperFigures');

files = betterdir(fullfile(resultDir));

fns = cell(length(files),1);

for f = 1:length(files)

    fns{f} = fullfile( ...
        files(f).folder, ...
        files(f).name, ...
        'NewEntropy_Raw.mat');

end

files = fns;


%% =========================================================
% Collect results
%% =========================================================

allData = [];


for fi = 1:length(files)

    fn = files{fi};

    S = load(fn);

    R = S.resultsAll;


    %% -----------------------------------------------------
    % Clear old effective dimensionality
    %% -----------------------------------------------------

    for ri = 1:length(R)

        R(ri).d_eff = [];

    end


    %% -----------------------------------------------------
    % Recalculate PI summary metrics
    %% -----------------------------------------------------

    for ri = 1:length(R)

        R50 = R(ri);


        % Restrict PI to desired history/horizon range

        R50.PI = R50.PI(1:10,1:3);

        R50.history = R50.history(1:10);


        PI = R50.PI;


        %% -------------------------------------------------
        % Raw PI maximum
        %% -------------------------------------------------

        R50.rawPI_max = ...
            max(PI(:),[],'omitnan');


        %% -------------------------------------------------
        % One-step PI
        %% -------------------------------------------------

        R50.rawPI_oneStep = ...
            PI(end,1);


        %% -------------------------------------------------
        % PI area
        %% -------------------------------------------------

        R50.rawPI_area = ...
            sum(PI(:),'omitnan');


        %% -------------------------------------------------
        % History timescale
        %% -------------------------------------------------

        PI_history = PI(:,1);

        if any(isfinite(PI_history))

            fracHistory = ...
                PI_history ./ ...
                max(PI_history,[],'omitnan');

            idx = find( ...
                fracHistory >= 0.95, ...
                1,'first');


            if isempty(idx)

                R50.historyTimescale_ms = NaN;

            else

                R50.historyTimescale_ms = ...
                    R50.history(idx)*1000;

            end

        else

            R50.historyTimescale_ms = NaN;

        end


        %% -------------------------------------------------
        % Future timescale
        %% -------------------------------------------------

        PI_future = PI(end,:);

        if any(isfinite(PI_future))

            fracFuture = ...
                PI_future ./ ...
                max(PI_future,[],'omitnan');

            idx = find( ...
                fracFuture <= 0.15, ...
                1,'first');


            if isempty(idx)

                R50.futureTimescale_ms = NaN;

            else

                R50.futureTimescale_ms = ...
                    R50.horizon(idx)*1000;

            end

        else

            R50.futureTimescale_ms = NaN;

        end


        %% -------------------------------------------------
        % Effective dimensionality
        %% -------------------------------------------------

        C = cumsum(R50.varExplained)/100;

        C = C(:)';

        p = diff([0 C]);


        if sum(p) > 0

            p = p ./ sum(p);

            R50.d_eff = ...
                1 ./ sum(p.^2);

        else

            R50.d_eff = NaN;

        end


        R(ri) = R50;

    end


    %% -----------------------------------------------------
    % Convert to long-form table
    %% -----------------------------------------------------

    for ri = 1:length(R)

        row.patient = ...
            string(R(ri).patient);

        row.cat = ...
            R(ri).cat;

        row.rawPI_max = ...
            R(ri).rawPI_max;

        row.pcaVariance = ...
            R(ri).pcaVariance;

        row.d_eff = ...
            R(ri).d_eff;

        row.rawPI_oneStep = ...
            R(ri).rawPI_oneStep;

        row.rawPI_area = ...
            R(ri).rawPI_area;

        row.normPI_max = ...
            R(ri).normPI_max;

        row.normPI_oneStep = ...
            R(ri).normPI_oneStep;

        row.historyTimescale_ms = ...
            R(ri).historyTimescale_ms;

        row.futureTimescale_ms = ...
            R(ri).futureTimescale_ms;

        row.nanFraction = ...
            R(ri).nanFraction;

        row.VarExp = {R(ri).varExplained};


        Trow = struct2table(row);

        allData = ...
            [allData; Trow];

    end

end


%% =========================================================
% Convert patient to categorical
%% =========================================================

allData.patient = ...
    categorical(allData.patient);


%% =========================================================
% Restrict PCA variance
%
% Same criterion used in previous analysis.
%
% IMPORTANT:
% Save the cutoff before exclusion so it can also be used
% when selecting the representative Pt XXX recordings
% for Tile 5.
%% =========================================================

pcaCutoff = ...
    prctile(allData.pcaVariance,90);


allData = ...
    allData(allData.pcaVariance <= pcaCutoff,:);


%% =========================================================
% MODEL DATASETS
%% =========================================================

validPI = ...
    ~ismissing(allData.cat) & ...
    isfinite(allData.rawPI_max);


validPI_PCA = ...
    ~ismissing(allData.cat) & ...
    isfinite(allData.rawPI_max) & ...
    isfinite(allData.pcaVariance);


validPCA = ...
    ~ismissing(allData.cat) & ...
    isfinite(allData.pcaVariance);


T_PI = ...
    allData(validPI,:);


T_PI_PCA = ...
    allData(validPI_PCA,:);


T_PCA = ...
    allData(validPCA,:);


%% =========================================================
% MODEL 1
%
% PImax ~ Severity
%% =========================================================

lme_PI = fitlme( ...
    T_PI, ...
    'rawPI_max ~ cat + (1|patient)');


%% =========================================================
% MODEL 2
%
% PImax ~ Severity + PCA variance
%% =========================================================

lme_PI_PCA = fitlme( ...
    T_PI_PCA, ...
    'rawPI_max ~ cat + pcaVariance + (1|patient)');


%% =========================================================
% MODEL 3
%
% PCA variance ~ Severity
%% =========================================================

lme_PCA = fitlme( ...
    T_PCA, ...
    'pcaVariance ~ cat + (1|patient)');


%% =========================================================
% Display models
%% =========================================================

disp('=================================================');
disp('MODEL 1: rawPI_max ~ cat + (1|patient)');
disp(lme_PI);

disp('=================================================');
disp('MODEL 2: rawPI_max ~ cat + pcaVariance + (1|patient)');
disp(lme_PI_PCA);

disp('=================================================');
disp('MODEL 3: pcaVariance ~ cat + (1|patient)');
disp(lme_PCA);


%% =========================================================
% STANDARDIZED EFFECTS
%
% LME models remain fit on RAW data.
%
% Standardization:
%
%   beta_std = beta_raw * SD(X)_within / SD(Y)_within
%
% The within-patient SDs are calculated separately for
% each model using the exact observations entering that model.
%% =========================================================


%% =========================================================
% MODEL 1
%
% PImax ~ Severity
%% =========================================================

idxCat1 = ...
    strcmp( ...
    lme_PI.CoefficientNames, ...
    'cat');


beta_PI_raw = ...
    lme_PI.Coefficients.Estimate(idxCat1);

p_PI = ...
    lme_PI.Coefficients.pValue(idxCat1);

CI_PI_raw = ...
    lme_PI.Coefficients.('Lower')(idxCat1);

CI_PI_raw_hi = ...
    lme_PI.Coefficients.('Upper')(idxCat1);


% ---------------------------------------------------------
% Within-patient SD of severity
% ---------------------------------------------------------

patients_PI = categories(T_PI.patient);

nPatients_PI = length(patients_PI);

Xsd = nan(nPatients_PI,1);
nX = nan(nPatients_PI,1);


for si = 1:nPatients_PI

    idx = ...
        T_PI.patient == patients_PI{si};

    x = ...
        T_PI.cat(idx);

    x = x(isfinite(x));

    if numel(x) >= 2

        Xsd(si) = std(x);
        nX(si) = numel(x);

    end

end


validSD = ...
    isfinite(Xsd) & ...
    nX > 1;


SDcat_PI = ...
    sqrt( ...
    sum( ...
    (nX(validSD)-1).*Xsd(validSD).^2) ./ ...
    sum(nX(validSD)-1));


% ---------------------------------------------------------
% Within-patient SD of PImax
% ---------------------------------------------------------

Ysd = nan(nPatients_PI,1);
nY = nan(nPatients_PI,1);


for si = 1:nPatients_PI

    idx = ...
        T_PI.patient == patients_PI{si};

    y = ...
        T_PI.rawPI_max(idx);

    y = y(isfinite(y));

    if numel(y) >= 2

        Ysd(si) = std(y);
        nY(si) = numel(y);

    end

end


validSD = ...
    isfinite(Ysd) & ...
    nY > 1;


SDPI_PI = ...
    sqrt( ...
    sum( ...
    (nY(validSD)-1).*Ysd(validSD).^2) ./ ...
    sum(nY(validSD)-1));


% ---------------------------------------------------------
% Standardized severity effect
% ---------------------------------------------------------

beta_PI_std = ...
    beta_PI_raw * ...
    SDcat_PI / SDPI_PI;


CI_PI_std = ...
    [ ...
    CI_PI_raw    * SDcat_PI / SDPI_PI, ...
    CI_PI_raw_hi * SDcat_PI / SDPI_PI];


%% =========================================================
% MODEL 2
%
% PImax ~ Severity + PCA variance
%% =========================================================

idxCat2 = ...
    strcmp( ...
    lme_PI_PCA.CoefficientNames, ...
    'cat');


idxPCA2 = ...
    strcmp( ...
    lme_PI_PCA.CoefficientNames, ...
    'pcaVariance');


% ---------------------------------------------------------
% Raw coefficients
% ---------------------------------------------------------

beta_PI_PCA_cat_raw = ...
    lme_PI_PCA.Coefficients.Estimate(idxCat2);

p_PI_PCA_cat = ...
    lme_PI_PCA.Coefficients.pValue(idxCat2);

CI_PI_PCA_cat = ...
    lme_PI_PCA.Coefficients.('Lower')(idxCat2);

CI_PI_PCA_cat_hi = ...
    lme_PI_PCA.Coefficients.('Upper')(idxCat2);


beta_PI_PCA_raw = ...
    lme_PI_PCA.Coefficients.Estimate(idxPCA2);

p_PI_PCA = ...
    lme_PI_PCA.Coefficients.pValue(idxPCA2);

CI_PI_PCA_raw = ...
    lme_PI_PCA.Coefficients.('Lower')(idxPCA2);

CI_PI_PCA_raw_hi = ...
    lme_PI_PCA.Coefficients.('Upper')(idxPCA2);


% ---------------------------------------------------------
% Within-patient SDs
%
% All SDs are calculated from T_PI_PCA because this is
% exactly the dataset used for Model 2.
% ---------------------------------------------------------

patients_PI_PCA = ...
    categories(T_PI_PCA.patient);

nPatients_PI_PCA = ...
    length(patients_PI_PCA);


% Severity SD

Xsd = nan(nPatients_PI_PCA,1);
nX = nan(nPatients_PI_PCA,1);


for si = 1:nPatients_PI_PCA

    idx = ...
        T_PI_PCA.patient == patients_PI_PCA{si};

    x = ...
        T_PI_PCA.cat(idx);

    x = x(isfinite(x));

    if numel(x) >= 2

        Xsd(si) = std(x);
        nX(si) = numel(x);

    end

end


validSD = ...
    isfinite(Xsd) & ...
    nX > 1;


SDcat_PI_PCA = ...
    sqrt( ...
    sum( ...
    (nX(validSD)-1).*Xsd(validSD).^2) ./ ...
    sum(nX(validSD)-1));


% PImax SD

Ysd = nan(nPatients_PI_PCA,1);
nY = nan(nPatients_PI_PCA,1);


for si = 1:nPatients_PI_PCA

    idx = ...
        T_PI_PCA.patient == patients_PI_PCA{si};

    y = ...
        T_PI_PCA.rawPI_max(idx);

    y = y(isfinite(y));

    if numel(y) >= 2

        Ysd(si) = std(y);
        nY(si) = numel(y);

    end

end


validSD = ...
    isfinite(Ysd) & ...
    nY > 1;


SDPI_PCA = ...
    sqrt( ...
    sum( ...
    (nY(validSD)-1).*Ysd(validSD).^2) ./ ...
    sum(nY(validSD)-1));


% PCA variance SD

PCA_sd = nan(nPatients_PI_PCA,1);
nPCA = nan(nPatients_PI_PCA,1);


for si = 1:nPatients_PI_PCA

    idx = ...
        T_PI_PCA.patient == patients_PI_PCA{si};

    y = ...
        T_PI_PCA.pcaVariance(idx);

    y = y(isfinite(y));

    if numel(y) >= 2

        PCA_sd(si) = std(y);
        nPCA(si) = numel(y);

    end

end


validSD = ...
    isfinite(PCA_sd) & ...
    nPCA > 1;


SDPCA_PI = ...
    sqrt( ...
    sum( ...
    (nPCA(validSD)-1).*PCA_sd(validSD).^2) ./ ...
    sum(nPCA(validSD)-1));


% ---------------------------------------------------------
% Standardized severity effect
% ---------------------------------------------------------

beta_PI_PCA_cat_std = ...
    beta_PI_PCA_cat_raw * ...
    SDcat_PI_PCA / SDPI_PCA;


CI_PI_PCA_cat_std = ...
    [ ...
    CI_PI_PCA_cat    * SDcat_PI_PCA / SDPI_PCA, ...
    CI_PI_PCA_cat_hi * SDcat_PI_PCA / SDPI_PCA];


% ---------------------------------------------------------
% Standardized PCA variance effect
% ---------------------------------------------------------

beta_PI_PCA_std = ...
    beta_PI_PCA_raw * ...
    SDPCA_PI / SDPI_PCA;


CI_PI_PCA_std = ...
    [ ...
    CI_PI_PCA_raw    * SDPCA_PI / SDPI_PCA, ...
    CI_PI_PCA_raw_hi * SDPCA_PI / SDPI_PCA];


%% =========================================================
% MODEL 3
%
% PCA variance ~ Severity
%% =========================================================

idxCat3 = ...
    strcmp( ...
    lme_PCA.CoefficientNames, ...
    'cat');


beta_PCA_raw = ...
    lme_PCA.Coefficients.Estimate(idxCat3);

p_PCA = ...
    lme_PCA.Coefficients.pValue(idxCat3);

CI_PCA_raw = ...
    lme_PCA.Coefficients.('Lower')(idxCat3);

CI_PCA_raw_hi = ...
    lme_PCA.Coefficients.('Upper')(idxCat3);


% ---------------------------------------------------------
% Within-patient SDs
%
% All from T_PCA because this is the dataset used for
% Model 3.
% ---------------------------------------------------------

patients_PCA = ...
    categories(T_PCA.patient);

nPatients_PCA = ...
    length(patients_PCA);


% Severity SD

Xsd = nan(nPatients_PCA,1);
nX = nan(nPatients_PCA,1);


for si = 1:nPatients_PCA

    idx = ...
        T_PCA.patient == patients_PCA{si};

    x = ...
        T_PCA.cat(idx);

    x = x(isfinite(x));

    if numel(x) >= 2

        Xsd(si) = std(x);
        nX(si) = numel(x);

    end

end


validSD = ...
    isfinite(Xsd) & ...
    nX > 1;


SDcat_PCA = ...
    sqrt( ...
    sum( ...
    (nX(validSD)-1).*Xsd(validSD).^2) ./ ...
    sum(nX(validSD)-1));


% PCA variance SD

PCA_sd = nan(nPatients_PCA,1);
nPCA = nan(nPatients_PCA,1);


for si = 1:nPatients_PCA

    idx = ...
        T_PCA.patient == patients_PCA{si};

    y = ...
        T_PCA.pcaVariance(idx);

    y = y(isfinite(y));

    if numel(y) >= 2

        PCA_sd(si) = std(y);
        nPCA(si) = numel(y);

    end

end


validSD = ...
    isfinite(PCA_sd) & ...
    nPCA > 1;


SDPCA_PCA = ...
    sqrt( ...
    sum( ...
    (nPCA(validSD)-1).*PCA_sd(validSD).^2) ./ ...
    sum(nPCA(validSD)-1));


% ---------------------------------------------------------
% Standardized severity effect
% ---------------------------------------------------------

beta_PCA_std = ...
    beta_PCA_raw * ...
    SDcat_PCA / SDPCA_PCA;


CI_PCA_std = ...
    [ ...
    CI_PCA_raw    * SDcat_PCA / SDPCA_PCA, ...
    CI_PCA_raw_hi * SDcat_PCA / SDPCA_PCA];


%% =========================================================
% ATTENUATION OF SEVERITY EFFECT
%
% Descriptive comparison, NOT a mediation estimate.
%% =========================================================

PI_severity_attenuation = ...
    100 * ...
    (beta_PI_std - beta_PI_PCA_cat_std) / ...
    beta_PI_std;


%% =========================================================
% Print standardized model results
%% =========================================================

fprintf('\n\n=================================================\n');
fprintf('STANDARDIZED EFFECTS\n');
fprintf('=================================================\n');

fprintf('\nModel 1: PI ~ Severity\n');

fprintf( ...
    '  beta = %.3f, p = %.3g\n', ...
    beta_PI_std, ...
    p_PI);


fprintf('\nModel 2: PI ~ Severity + PCA variance\n');

fprintf( ...
    '  Severity beta = %.3f, p = %.3g\n', ...
    beta_PI_PCA_cat_std, ...
    p_PI_PCA_cat);

fprintf( ...
    '  PCA variance beta = %.3f, p = %.3g\n', ...
    beta_PI_PCA_std, ...
    p_PI_PCA);


fprintf('\nModel 3: PCA variance ~ Severity\n');

fprintf( ...
    '  beta = %.3f, p = %.3g\n', ...
    beta_PCA_std, ...
    p_PCA);


fprintf('\nSeverity beta attenuation after PCA adjustment:\n');

fprintf( ...
    '  %.1f%%\n', ...
    PI_severity_attenuation);


%% =========================================================
% PLOTTING SETTINGS
%% =========================================================

ptSize = 40;
lnWidth = 1;
fontSize = 12;
titleFontSize = 16;

presentColor = [213 94 0] / 255;
pastColor    = [0 158 115] / 255;

nPatients = length(allpts);

cmap = lines(nPatients);


%% =========================================================
% CREATE FIGURE
%
% 4 x 4 layout
%
% Row 1:
%   [ Tile 1: BLANK ] [ Tile 4: PCA ]
%
% Row 2:
%   [ Tile 2: SIGNAL ] [ Tile 4: PCA ]
%
% Row 3:
%   [ Tile 3: PI ] [ Tile 5: PC SPECTRUM ]
%
% Row 4:
%   [ Tile 3: PI ] [ Tile 6: FOREST ]
%
%
% Exact spans:
%
% Tile 1 = row 1, columns 1:2
% Tile 2 = row 2, columns 1:2
% Tile 3 = rows 3:4, columns 1:2
% Tile 4 = rows 1:2, columns 3:4
% Tile 5 = row 3, columns 3:4
% Tile 6 = row 4, columns 3:4
%% =========================================================

fig = figure( ...
    'Color','w', ...
    'Position',[100 50 1200 1000]);


tl = tiledlayout( ...
    fig, ...
    4,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');





%% =========================================================
% TILE 2
%
% REPRESENTATIVE NEURAL TRACES
%
% Row 2, columns 1:2.
%% =========================================================

ax = nexttile(tl,2,[2 1]);

hold(ax,'on');

%% --------------------------------------------------------- % Load representative signal %% --------------------------------------------------------- 
fn = ''; %example file
reginfo = load(betterLoadSave([USER_DIR,'/Info/',examplePatient,'_Regions.mat']));
keepreg = [1,3,4,5,7,9,10,11]; 
 
regnames = reginfo.region_names;
        for rq = 1:length(regnames)
            regnames{rq} = char(strrep(regnames{rq},"'",""));
        end
        indskp = find(ismember(reginfo.region_labels,keepreg));
        regnames = regnames(indskp);

Ssignal = load(fn); 
if ~isfield(Ssignal,'data') 
    error( ... 
        'Representative signal file does not contain variable "data".'); 
end 
data = Ssignal.data; 



            datachans = Ssignal.chanNames(2:end)';
            keepelecs = find(ismember(datachans,regnames));
data = data(:,keepelecs);
if isfield(Ssignal,'timevec') 
    timevec = Ssignal.timevec(:); 
elseif isfield(Ssignal,'sfreq') 
    
    timevec = ... 
        (0:size(data,1)-1)' ./ Ssignal.sfreq; 
else 
    error( ... 
        'Representative signal file contains neither "timevec" nor "sfreq".'); 
end 
% Ensure time is the first dimension 
if size(data,1) ~= length(timevec) && ... 
        size(data,2) == length(timevec) 
    data = data'; 
end
%% ---------------------------------------------------------
% PCA of representative signal
%% ---------------------------------------------------------

% Ensure time is the first dimension
if size(data,1) ~= length(timevec) && ...
        size(data,2) == length(timevec)

    data = data';

end

% Z-score each channel across time
data = zscore(data,0,1);

% PCA across channels
[~,score,~,~,explained] = pca(data);

% Retain top 10 PCs
nPC = min(5,size(score,2));

score = score(:,1:nPC);


%% ---------------------------------------------------------
% Plot top 10 PCs
%% ---------------------------------------------------------

for j = 1:nPC

    ts = 2*rescale(score(:,j));

    shift = j - 0.5;

    plot( ...
        ax, ...
        timevec, ...
        ts + shift, ...
        'k', ...
        'LineWidth',2);

end


%% ---------------------------------------------------------
% Get plot limits
%% ---------------------------------------------------------

xlim(ax,[29.7 30.1]);

yl = ylim(ax);

ypatch = ...
    [yl(1) yl(1) yl(2) yl(2)];

ylim([1 nPC+1.5])


%% ---------------------------------------------------------
% Current window: 30 to 30.005 s
%% ---------------------------------------------------------

xpatch = ...
    [30 30.005 30.005 30];


hp = patch( ...
    ax, ...
    xpatch, ...
    ypatch, ...
    presentColor, ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none');


uistack(hp,'bottom');


%% ---------------------------------------------------------
% Past windows: 5 ms each, back to 200 ms
%% ---------------------------------------------------------

nPatches = 40;

startAlpha = 0.5;

endAlpha = 0.03;


for k = 1:nPatches

    tEnd = ...
        30 - (k-1)*0.005;

    tStart = ...
        30 - k*0.005;


    xpatch = ...
        [tStart tEnd tEnd tStart];


    alphaVal = ...
        startAlpha - ...
        (k-1) * ...
        (startAlpha-endAlpha) / ...
        (nPatches-1);


    hp = patch( ...
        ax, ...
        xpatch, ...
        ypatch, ...
        pastColor, ...
        'FaceAlpha',alphaVal, ...
        'EdgeColor','none');


    uistack(hp,'bottom');

end


%% ---------------------------------------------------------
% Formatting
%% ---------------------------------------------------------

xlim(ax,[29.7 30.1]);

set( ...
    ax, ...
    'XTick',[], ...
    'YTick',[], ...
    'FontSize',fontSize, ...
    'Box','off');


ylabel( ...
    ax, ...
    'Principal Components', ...
    'FontSize',fontSize);


xlabel( ...
    ax, ...
    'Time', ...
    'FontSize',fontSize);

%% =========================================================
% TILE 1
%
% BLANK
%
% Upper-left.
%% =========================================================

nexttile(tl,1,[2 1]);
plot(data(:,end),'k');
ylim([-25,25]);
axis off
box off


fig2 = figure( ...
    'Color','w', ...
    'Position',[100 50 1200 1000]);


tl2 = tiledlayout( ...
    fig2, ...
    4,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');
nexttile(tl2,1,[2 1]);
imagesc(data');
xlabel('Time')
ylabel('Bipolar Channel Voltage (Z-score)')
xticklabels([]);
yticklabels([]);
clim([-4 4]);


%% =========================================================
% TILE 3
%
% RAW PREDICTIVE INFORMATION
%
% Rows 3:4, columns 1:2.
%
% IMPORTANT:
% Visualization uses WITHIN-SUBJECT Z-SCORES exactly as
% in the original code.
%
% The LME statistics remain based on RAW data.
%% =========================================================

ax = nexttile(tl,9,[2 2]);

hold(ax,'on');


subjects = categories(T_PI.patient);


XzAll = [];

YzAll = [];


% Store subject handles for optional legend
hSubjects = gobjects(nPatients,1);


for si = 1:nPatients

    idx = ...
        T_PI.patient == subjects{si};


    if nnz(idx) < 2

        continue

    end


    xRaw = ...
        T_PI.cat(idx);

    yRaw = ...
        T_PI.rawPI_max(idx);


    %% -----------------------------------------------------
    % Within-subject z-score
    %% -----------------------------------------------------

    muX = mean(xRaw);

    sdX = std(xRaw);


    muY = mean(yRaw);

    sdY = std(yRaw);


    if sdX == 0 || ...
            sdY == 0 || ...
            ~isfinite(sdX) || ...
            ~isfinite(sdY)

        continue

    end


    x = ...
        (xRaw-muX) ./ sdX;

    y = ...
        (yRaw-muY) ./ sdY;


    %% -----------------------------------------------------
    % Add to pooled z-scored dataset
    %% -----------------------------------------------------

    XzAll = ...
        [XzAll; x(:)];

    YzAll = ...
        [YzAll; y(:)];


    %% -----------------------------------------------------
    % Scatter
    %% -----------------------------------------------------

    hSubjects(si) = scatter( ...
        ax, ...
        x, ...
        y, ...
        ptSize, ...
        'filled', ...
        'MarkerFaceColor',cmap(si,:), ...
        'MarkerFaceAlpha',0.5, ...
        'MarkerEdgeColor','none');


    %% -----------------------------------------------------
    % Subject-specific regression
    %% -----------------------------------------------------

    p = ...
        polyfit(x,y,1);


    xfit = ...
        linspace(min(x),max(x),100);


    yfit = ...
        polyval(p,xfit);


    plot( ...
        ax, ...
        xfit, ...
        yfit, ...
        '--', ...
        'Color',cmap(si,:), ...
        'LineWidth',1.5);

end


%% ---------------------------------------------------------
% Pooled visualization regression
%
% Fit to the concatenated within-subject z-scored values.
%% ---------------------------------------------------------

mdlPlot_PI = ...
    fitlm(XzAll,YzAll);


xfitPlot = ...
    linspace( ...
    min(XzAll), ...
    max(XzAll), ...
    100)';


[yfitPlot,yCI] = ...
    predict( ...
    mdlPlot_PI, ...
    xfitPlot);


rPlot_PI = ...
    corr( ...
    XzAll, ...
    YzAll, ...
    'Rows','complete');


fill( ...
    ax, ...
    [xfitPlot;flipud(xfitPlot)], ...
    [yCI(:,1);flipud(yCI(:,2))], ...
    [0.8 0.8 0.8], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.3);


plot( ...
    ax, ...
    xfitPlot, ...
    yfitPlot, ...
    'k-', ...
    'LineWidth',2);


%% ---------------------------------------------------------
% PI annotation
%% ---------------------------------------------------------

text( ...
    ax, ...
    0.98, ...
    0.98, ...
    { ...
    sprintf( ...
    'Model 1: β = %.3f, p = %.3g', ...
    beta_PI_std, ...
    p_PI), ...
    sprintf( ...
    'Model 2: β_{severity} = %.3f, p = %.3g', ...
    beta_PI_PCA_cat_std, ...
    p_PI_PCA_cat), ...
    sprintf( ...
    'Pooled r = %.3f', ...
    rPlot_PI)}, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);


%% ---------------------------------------------------------
% Labels
%% ---------------------------------------------------------

xlabel( ...
    ax, ...
    {'Depression Severity','(Z Score)'}, ...
    'FontSize',fontSize);


ylabel( ...
    ax, ...
    {'Predictive Information','(Z Score)'}, ...
    'FontSize',fontSize);


title( ...
    ax, ...
    'Raw Predictive Information', ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal');


xlim(ax,[-2.5 2.5]);

axis(ax,'square');

box(ax,'off');


%% =========================================================
% TILE 4
%
% PCA VARIANCE EXPLAINED
%
% Rows 1:2, columns 3:4.
%
% Visualization uses WITHIN-SUBJECT Z-SCORES.
%% =========================================================

ax = nexttile(tl,3,[2 2]);

hold(ax,'on');


subjects = categories(T_PCA.patient);


XzAll = [];

YzAll = [];


hSubjectsPCA = gobjects(nPatients,1);


for si = 1:nPatients

    idx = ...
        T_PCA.patient == subjects{si};


    if nnz(idx) < 2

        continue

    end


    xRaw = ...
        T_PCA.cat(idx);

    yRaw = ...
        T_PCA.pcaVariance(idx);


    %% -----------------------------------------------------
    % Within-subject z-score
    %% -----------------------------------------------------

    muX = mean(xRaw);

    sdX = std(xRaw);


    muY = mean(yRaw);

    sdY = std(yRaw);


    if sdX == 0 || ...
            sdY == 0 || ...
            ~isfinite(sdX) || ...
            ~isfinite(sdY)

        continue

    end


    x = ...
        (xRaw-muX) ./ sdX;

    y = ...
        (yRaw-muY) ./ sdY;


    %% -----------------------------------------------------
    % Add to pooled z-scored dataset
    %% -----------------------------------------------------

    XzAll = ...
        [XzAll; x(:)];

    YzAll = ...
        [YzAll; y(:)];


    %% -----------------------------------------------------
    % Scatter
    %% -----------------------------------------------------

    hSubjectsPCA(si) = scatter( ...
        ax, ...
        x, ...
        y, ...
        ptSize, ...
        'filled', ...
        'MarkerFaceColor',cmap(si,:), ...
        'MarkerFaceAlpha',0.5, ...
        'MarkerEdgeColor','none');


    %% -----------------------------------------------------
    % Subject-specific regression
    %% -----------------------------------------------------

    p = ...
        polyfit(x,y,1);


    xfit = ...
        linspace(min(x),max(x),100);


    yfit = ...
        polyval(p,xfit);


    plot( ...
        ax, ...
        xfit, ...
        yfit, ...
        '--', ...
        'Color',cmap(si,:), ...
        'LineWidth',1.5);

end


%% ---------------------------------------------------------
% Pooled visualization regression
%% ---------------------------------------------------------

mdlPlot_PCA = ...
    fitlm(XzAll,YzAll);


xfitPlot = ...
    linspace( ...
    min(XzAll), ...
    max(XzAll), ...
    100)';


[yfitPlot,yCI] = ...
    predict( ...
    mdlPlot_PCA, ...
    xfitPlot);


rPlot_PCA = ...
    corr( ...
    XzAll, ...
    YzAll, ...
    'Rows','complete');


fill( ...
    ax, ...
    [xfitPlot;flipud(xfitPlot)], ...
    [yCI(:,1);flipud(yCI(:,2))], ...
    [0.8 0.8 0.8], ...
    'EdgeColor','none', ...
    'FaceAlpha',0.3);


plot( ...
    ax, ...
    xfitPlot, ...
    yfitPlot, ...
    'k-', ...
    'LineWidth',2);


%% ---------------------------------------------------------
% PCA annotation
%% ---------------------------------------------------------

text( ...
    ax, ...
    0.98, ...
    0.98, ...
    { ...
    sprintf( ...
    'Model 3: β = %.3f, p = %.3g', ...
    beta_PCA_std, ...
    p_PCA), ...
    sprintf( ...
    'Pooled r = %.3f', ...
    rPlot_PCA)}, ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',fontSize);


%% ---------------------------------------------------------
% Labels
%% ---------------------------------------------------------

xlabel( ...
    ax, ...
    {'Depression Severity','(Z Score)'}, ...
    'FontSize',fontSize);


ylabel( ...
    ax, ...
    {'PCA Variance Explained','(Z Score)'}, ...
    'FontSize',fontSize);


title( ...
    ax, ...
    'PCA Variance Explained', ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal');


xlim(ax,[-2.5 2.5]);

axis(ax,'square');

box(ax,'off');


%% ---------------------------------------------------------
% Subject legend
%% ---------------------------------------------------------

validHandles = ...
    isgraphics(hSubjectsPCA);


legend( ...
    ax, ...
    hSubjectsPCA(validHandles), ...
    alllabels(validHandles), ...
    'Location','southeast', ...
    'Box','off');


%% =========================================================
% TILE 5
%
% TOP 20 PRINCIPAL COMPONENTS
%
% Least vs most depressed Pt XXX recording.
%
% Row 3, columns 3:4.
%
% IMPORTANT:
% The same global 90th-percentile PCA exclusion used above
% is applied BEFORE identifying the least/most depressed
% recordings.
%% =========================================================

ax = nexttile(tl,11,[1 2]);

hold(ax,'on');


%% =========================================================
% PCA VARIANCE SPECTRA ACROSS SUBJECTS
%
% For each subject:
%   - identify least-depressed surviving recording
%   - identify most-depressed surviving recording
%
% Then:
%   - calculate mean variance explained for each PC
%   - calculate SEM across subjects
%   - plot least vs most depressed recordings
%% =========================================================

%% ---------------------------------------------------------
% Parameters
%% ---------------------------------------------------------

nPC = 20;

fontSize = 9;
titleFontSize = 12;
lnWidth = 1;

%% ---------------------------------------------------------
% Subjects
%% ---------------------------------------------------------

subjects = categories(allData.patient);
nSubjects = length(subjects);

%% ---------------------------------------------------------
% Preallocate
%
% Rows = subjects
% Columns = PCs
%% ---------------------------------------------------------

varExpLeast = nan(nSubjects,nPC);
varExpMost  = nan(nSubjects,nPC);

catLeast = nan(nSubjects,1);
catMost  = nan(nSubjects,1);

%% ---------------------------------------------------------
% Select least- and most-depressed recording within
% each subject
%% ---------------------------------------------------------

for si = 1:nSubjects

    thisSubject = subjects{si};

    idxSubject = ...
        allData.patient == thisSubject;

    Tsub = ...
        allData(idxSubject,:);

    % Remove rows without CAT-DI
    valid = ...
        isfinite(Tsub.cat);

    Tsub = ...
        Tsub(valid,:);

    if isempty(Tsub)

        continue

    end

    %% -----------------------------------------------------
    % Least depressed recording
    %% -----------------------------------------------------

    [catLeast(si),idxLeast] = ...
        min(Tsub.cat);

    rowLeast = ...
        Tsub(idxLeast,:);

    ve = ...
        rowLeast.VarExp{1};

    ve = ...
        ve(:);

    nUse = ...
        min(nPC,length(ve));

    varExpLeast(si,1:nUse) = ...
        ve(1:nUse);


    %% -----------------------------------------------------
    % Most depressed recording
    %% -----------------------------------------------------

    [catMost(si),idxMost] = ...
        max(Tsub.cat);

    rowMost = ...
        Tsub(idxMost,:);

    ve = ...
        rowMost.VarExp{1};

    ve = ...
        ve(:);

    nUse = ...
        min(nPC,length(ve));

    varExpMost(si,1:nUse) = ...
        ve(1:nUse);

end

%% ---------------------------------------------------------
% Mean and SEM across subjects
%% ---------------------------------------------------------

meanLeast = ...
    mean(varExpLeast,1,'omitnan');

meanMost = ...
    mean(varExpMost,1,'omitnan');

nLeast = ...
    sum(isfinite(varExpLeast),1);

nMost = ...
    sum(isfinite(varExpMost),1);

semLeast = ...
    std(varExpLeast,0,1,'omitnan') ./ sqrt(nLeast);

semMost = ...
    std(varExpMost,0,1,'omitnan') ./ sqrt(nMost);

%% ---------------------------------------------------------
% Plot
%% ---------------------------------------------------------

blue = [0.230 0.299 0.754];
red  = [0.706 0.016 0.150];

% Dimmed versions for PCs 11–20
blueDim = 0.55 * blue + 0.45 * [1 1 1];
redDim  = 0.55 * red  + 0.45 * [1 1 1];

pc = 1:nPC;

b = bar( ...
    ax, ...
    pc, ...
    [meanLeast(:) meanMost(:)], ...
    'grouped');

hold(ax,'on');

%% ---------------------------------------------------------
% Bar colors
%% ---------------------------------------------------------

% PCs 1–10
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';

b(1).CData = ...
    repmat(red,nPC,1);

b(2).CData = ...
    repmat(blue,nPC,1);

% PCs 11–20
if nPC > 10

    b(1).CData(11:end,:) = ...
        repmat(redDim,nPC-10,1);

    b(2).CData(11:end,:) = ...
        repmat(blueDim,nPC-10,1);

end

%% ---------------------------------------------------------
% Error bars
%% ---------------------------------------------------------

xLeast = ...
    b(1).XEndPoints;

xMost = ...
    b(2).XEndPoints;

errorbar( ...
    xLeast, ...
    meanLeast, ...
    semLeast, ...
    'k', ...
    'LineStyle','none', ...
    'LineWidth',lnWidth);

errorbar( ...
    xMost, ...
    meanMost, ...
    semMost, ...
    'k', ...
    'LineStyle','none', ...
    'LineWidth',lnWidth);

%% ---------------------------------------------------------
% Separation between PC 10 and PC 11
%% ---------------------------------------------------------

xline( ...
    ax, ...
    10.5, ...
    '--k', ...
    'LineWidth',1);

%% ---------------------------------------------------------
% Formatting
%% ---------------------------------------------------------

xlabel( ...
    ax, ...
    'Principal Component', ...
    'FontSize',fontSize);

ylabel( ...
    ax, ...
    'Variance Explained (%)', ...
    'FontSize',fontSize);

title( ...
    ax, ...
    'PCA variance spectra across subjects', ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal');

legend( ...
    ax, ...
    { ...
    'Least depressed recording', ...
    'Most depressed recording'}, ...
    'Location','northeast', ...
    'Box','off');

set( ...
    ax, ...
    'FontSize',fontSize, ...
    'Box','off', ...
    'TickDir','out');

xticks(ax,1:nPC);

xlim(ax,[0.5 nPC+0.5]);

hold(ax,'off');
%% =========================================================
% TILE 6
%
% FOREST PLOT
%
% Rows:
%
%   Model 1:
%       PI ~ Severity
%
%   Model 2:
%       PI ~ Severity + PCA
%       Severity coefficient
%
%   Model 3:
%       PCA ~ Severity
%
% Row 4, columns 3:4.
%% =========================================================

ax = nexttile(tl,15,[1 2]);

hold(ax,'on');


forestBeta = [ ...
    beta_PI_std; ...
    beta_PI_PCA_cat_std; ...
    beta_PCA_std];


forestCIlo = [ ...
    CI_PI_std(1); ...
    CI_PI_PCA_cat_std(1); ...
    CI_PCA_std(1)];


forestCIhi = [ ...
    CI_PI_std(2); ...
    CI_PI_PCA_cat_std(2); ...
    CI_PCA_std(2)];


forestP = [ ...
    p_PI; ...
    p_PI_PCA_cat; ...
    p_PCA];


forestLabels = { ...
    'Model 1'; ...
    'Model 2'; ...
    'Model 3'};


ypos = ...
    1:length(forestBeta);


%% ---------------------------------------------------------
% Zero line
%% ---------------------------------------------------------

plot( ...
    ax, ...
    [0 0], ...
    [0.5 length(ypos)+0.5], ...
    'k--', ...
    'LineWidth',1);


%% ---------------------------------------------------------
% Confidence intervals
%% ---------------------------------------------------------

for i = 1:length(ypos)

    if forestP(i) < 0.05

        lineColor = [0 0 0];

    else

        lineColor = [0.6 0.6 0.6];

    end


    plot( ...
        ax, ...
        [forestCIlo(i) forestCIhi(i)], ...
        [ypos(i) ypos(i)], ...
        '-', ...
        'Color',lineColor, ...
        'LineWidth',2);


    plot( ...
        ax, ...
        forestBeta(i), ...
        ypos(i), ...
        'o', ...
        'MarkerFaceColor',lineColor, ...
        'MarkerEdgeColor',lineColor, ...
        'MarkerSize',6);

end


%% ---------------------------------------------------------
% Labels
%% ---------------------------------------------------------

set( ...
    ax, ...
    'YTick',ypos, ...
    'YTickLabel',forestLabels, ...
    'YDir','reverse');


xlabel( ...
    ax, ...
    'Standardized LME Effect of Severity', ...
    'FontSize',fontSize);


title( ...
    ax, ...
    'LME effect estimates', ...
    'FontSize',titleFontSize, ...
    'FontWeight','normal');


box(ax,'off');


%% ---------------------------------------------------------
% Add beta and p-value text
%% ---------------------------------------------------------

xLimits = ...
    xlim(ax);


xText = ...
    xLimits(2);


for i = 1:length(ypos)

    if i == 2

        labelText = sprintf( ...
            'β_{severity} = %.3f, p = %.3g', ...
            forestBeta(i), ...
            forestP(i));

    else

        labelText = sprintf( ...
            'β = %.3f, p = %.3g', ...
            forestBeta(i), ...
            forestP(i));

    end


    text( ...
        ax, ...
        xText, ...
        ypos(i), ...
        labelText, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10);

end


%% =========================================================
% GENERAL FIGURE FORMATTING
%% =========================================================

set( ...
    findall(fig,'Type','axes'), ...
    'FontName','Arial');


set( ...
    findall(fig,'Type','text'), ...
    'FontName','Arial');


%% =========================================================
% SAVE
%% =========================================================
statsDir = [USER_DIR,'/StatTables'];

exportFigure6Statistics(T_PI,T_PI_PCA,T_PCA,allData,pcaCutoff,alllabels, ...
    fullfile(statsDir,'Figure6_Statistics.xlsx'));




exportgraphics( ...
    fig, ...
    betterLoadSave(fullfile(outDir,'Fig6_Predictability.svg')), ...
    'ContentType','vector');

exportgraphics( ...
    fig2, ...
    betterLoadSave(fullfile(outDir,'Fig6_Predictability_raster.png')));

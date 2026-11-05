function results = predictiveInformationHorizon(data, fs, targetFs, maxHistory, maxHorizon, nPC)
% PREDICTIVEINFORMATIONHORIZON
%
% Estimate Gaussian predictive information landscape:
%
% I(history;future)
%
% INPUTS:
%   data        : nElec x nTime
%   fs          : original sampling frequency
%   targetFs    : resampled frequency
%   maxHistory  : maximum history length (samples)
%   maxHorizon  : maximum prediction horizon (samples)
%   nPC         : number of PCA dimensions
%
% OUTPUT:
%   results.PI                  : predictive information matrix
%                                history x horizon
%
%   Summary metrics:
%   results.rawPI_max
%   results.rawPI_oneStep
%   results.rawPI_area
%   results.normPI_max
%   results.normPI_oneStep
%   results.historyTimescale_ms
%   results.futureTimescale_ms
%
%   results.pcaVariance
%   results.nanFraction


%% Downsample

if targetFs < fs
    data = resample(data',targetFs,fs)';
    fs = targetFs;
end


%% Normalize channels

data = zscore(data,0,2);


%% PCA

[~,score,~,~,explained] = pca(data');

if nPC > size(score,2)
    nPC = size(score,2);
end

score = score(:,1:nPC);

results.pcaVariance = single(sum(explained(1:nPC)));


%% Compute predictive information matrix

PI = nan(maxHistory,maxHorizon);

nTime = size(score,1);


for k = 1:maxHistory

    for h = 1:maxHorizon


        % Valid time indices
        idx = (k+1):(nTime-h);


        % Future state

        future = score(idx+h,:);


        % Past history

        past = [];

        for lag = k:-1:1
            past = [past score(idx-lag+1,:)];
        end


        % Remove invalid samples

        valid = all(isfinite(past),2) & ...
                all(isfinite(future),2);

        past = past(valid,:);
        future = future(valid,:);


        if size(past,1) < size(past,2)+5
            continue
        end


        %% Covariance matrices

        Cpast = cov(past);
        Cfuture = cov(future);
        Cjoint = cov([past future]);


        % Scale-dependent regularization

        Cpast = regularizeCov(Cpast);
        Cfuture = regularizeCov(Cfuture);
        Cjoint = regularizeCov(Cjoint);


        %% Gaussian mutual information

        logDetPast   = logdet(Cpast);
        logDetFuture = logdet(Cfuture);
        logDetJoint  = logdet(Cjoint);


        if all(isfinite([logDetPast logDetFuture logDetJoint]))

            PI(k,h) = 0.5/log(2) * ...
                (logDetPast + logDetFuture - logDetJoint);

        end

    end

end


%% Store primary outputs

results.PI = single(PI);

results.history = single((1:maxHistory)/fs);

results.horizon = single((1:maxHorizon)/fs);

results.fs = single(fs);

results.varExplained = single(explained);
%% Quality metric

results.nanFraction = single(mean(isnan(PI(:))));


%% Summary predictive information metrics


% Raw PI

results.rawPI_max = single(max(PI(:),[],'omitnan'));

results.rawPI_oneStep = single(PI(end,1));

results.rawPI_area = single(sum(PI(:),'omitnan'));



%% Normalize predictive information

% Entropy of future state

future = score(2:end,:);

Cfuture = cov(future);

Cfuture = regularizeCov(Cfuture);


logDetFuture = logdet(Cfuture);


if isfinite(logDetFuture)

    Hfuture = 0.5/log(2) * ...
        (nPC*log(2*pi*exp(1)) + logDetFuture);

    results.normPI_max = ...
        single(results.rawPI_max/Hfuture);

    results.normPI_oneStep = ...
        single(results.rawPI_oneStep/Hfuture);

else

    results.normPI_max = NaN;
    results.normPI_oneStep = NaN;

end



%% Predictive timescales


% History saturation:
% history required for 95% of one-step PI

PI_history = PI(:,1);

if any(isfinite(PI_history))

    fracHistory = PI_history ./ max(PI_history,[],'omitnan');

    idx = find(fracHistory >= 0.95,1,'first');

    if isempty(idx)
        results.historyTimescale_ms = NaN;
    else
        results.historyTimescale_ms = ...
            single(results.history(idx)*1000);
    end

else

    results.historyTimescale_ms = NaN;

end



% Future prediction horizon:
% horizon where PI falls below 5%

PI_future = PI(end,:);

if any(isfinite(PI_future))

    fracFuture = PI_future ./ max(PI_future,[],'omitnan');

    idx = find(fracFuture <= 0.05,1,'first');

    if isempty(idx)

        results.futureTimescale_ms = NaN;

    else

        results.futureTimescale_ms = ...
            single(results.horizon(idx)*1000);

    end

else

    results.futureTimescale_ms = NaN;

end


end

function result = pfunc_04_Inferrer_ML_Group(binarizedData, h_star,J_star, lr_rate, iter_max, stopping)
% Joint MLE of shared h*, J*, and per-recording beta_r with hard mask
% binarizedData: cell array {1 x R}, each cell = N x T_r (binary +1/-1)
% Outputs:
%   result.h_star: N x 1 shared h
%   result.J_star: N x N shared J
%   result.beta: R x 1 effective temperature per recording

if nargin < 4, lr_rate = 0.01; end
if nargin < 5, iter_max = 1e5; end
if nargin < 6, stopping = 5e-6; end
lr_rate_beta = lr_rate/5;
R = numel(binarizedData);
N = size(binarizedData{1},1);

% Initialize shared structure
%h_star = zeros(N,1);
%J_star = zeros(N);

% Initialize beta
beta = ones(R,1);

% Hard mask initialization
h_mask = true(N,1);
J_mask = true(N);

% Warm-up iteration for zeroing
iter_warmup = 1000;
threshold_factor = 0.1; % fraction of typical J* scale

% Precompute state list


% stateList = mfunc_VectorList(N);

for iter = 1:iter_max
    dh_total = zeros(N,1);
    dJ_total = zeros(N);

    % -----------------------------
    % Recompute data mean and correlation for current h*, J* (optional)
    % -----------------------------
    dataMean = zeros(N,R);
    dataCorr = zeros(N,N,R);
    dataEnergyMean = zeros(R,1);
    for r = 1:R
        X = binarizedData{r};
        T = size(X,2);
        dataMean(:,r) = mean(X,2);
        dataCorr(:,:,r) = (X*X')/T;

        % Compute data mean energy under current h*, J*
        E_data = zeros(1,T);
        for t = 1:T
            sigma = X(:,t);
            E_data(t) = -0.5 * sigma' * J_star * sigma - h_star' * sigma;
        end
        dataEnergyMean(r) = mean(E_data);
    end

    % -----------------------------
    % Update shared h* and J*
    % -----------------------------
    for r = 1:R
        [modelMean, modelCorr] = mfunc_ModelMeanCorrelation(h_star, J_star);

        % Weighted gradient for h
        dh_total = dh_total + beta(r)*(dataMean(:,r)-modelMean);

        % Weighted gradient for J, zero diagonal
        dJ_r = beta(r)*(dataCorr(:,:,r)-modelCorr);
        dJ_r = dJ_r - diag(diag(dJ_r));
        dJ_total = dJ_total + dJ_r;
    end

    % Update shared parameters with hard mask
    h_star = h_star + lr_rate*(dh_total .* h_mask);
    J_star = J_star + lr_rate*(dJ_total .* J_mask);

    % Apply zeroing and set mask after warm-up
    if iter == iter_warmup
        J_ut = triu(abs(J_star),1);
        J_scale = median(J_ut(:));
        thresh = threshold_factor * J_scale;

        h_mask = abs(h_star) >= thresh;
        J_mask = abs(J_star) >= thresh;

        % Apply zeroing
        h_star(~h_mask) = 0;
        J_star(~J_mask) = 0;
    end

    % -----------------------------
    % -----------------------------
    % Update beta per recording
    % -----------------------------
    for r = 1:R
        % Scale shared parameters by beta(r) for probability computation
        h_eff = beta(r) * h_star;
        J_eff = beta(r) * J_star;

        % Compute Boltzmann probabilities under current beta
        probModel_r = mfunc_StateProb(h_eff, J_eff);

        % Compute expected energy under original h*, J*
        modelEnergies = mfunc_Energy(h_star, J_star);
        E_modelMean_r = sum(probModel_r .* modelEnergies);

        % Gradient and update
        grad_beta = E_modelMean_r - dataEnergyMean(r);
        beta(r) = max(beta(r) + lr_rate_beta * grad_beta, 1e-8); % enforce positivity
    end


    % -----------------------------
    % Check convergence
    % -----------------------------
    change = sqrt(norm(dJ_total,'fro')^2 + norm(dh_total)^2)/N/(N+1);
    if change < stopping
        fprintf('Converged in %d iterations\n', iter);
        break;
    end
end

% -----------------------------
% Return results
% -----------------------------
result.h_star = h_star;
result.J_star = J_star;
result.beta = beta;
end

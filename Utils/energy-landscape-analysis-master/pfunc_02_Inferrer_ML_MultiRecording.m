function result = pfunc_02_Inferrer_ML_MultiRecording(dataCell,ablation, lr_rate, iter_max, stopping, eta_beta, lambda_h, lambda_J, check_interval)
%==========================================================================
% Multi-recording Ising MLE with per-recording beta updates
% Exact log-likelihood per recording used as convergence criterion
%==========================================================================
% 2025/10/19+ by S. V. Rajesh
%
% INPUTS:
%   dataCell       : cell array of R recordings; each [N x T_r] with -1/+1
%   lr_rate        : learning rate for h,J updates (default 0.05)
%   iter_max       : max iterations (default 1e5)
%   stopping       : relative log-likelihood improvement threshold (default 1e-6)
%   eta_beta       : learning rate for beta updates (default 0.01)
%   lambda_h       : L2 shrinkage for small h (default 1e-3)
%   check_interval : iterations between log-likelihood checks (default 50)
%
% OUTPUT:
%   result.h          -> inferred local fields [N x 1]
%   result.J          -> inferred couplings [N x N]
%   result.beta       -> inferred inverse temperatures [R x 1]
%   result.logL_hist  -> log-likelihood trace
%   result.iter       -> final iteration count
%==========================================================================
if nargin < 2 || isempty(ablation), ablation = 0; end
if nargin < 3 || isempty(lr_rate), lr_rate = 0.05; end
if nargin < 4 || isempty(iter_max), iter_max = 1e5; end
if nargin < 5 || isempty(stopping), stopping = 1e-6; end
if nargin < 6 || isempty(eta_beta), eta_beta = 0.01; end
if nargin < 7 || isempty(lambda_h), lambda_h = 1e-2; end
if nargin < 8 || isempty(lambda_J), lambda_J = 1e-2; end
if nargin < 9 || isempty(check_interval), check_interval = 50; end

R = numel(dataCell);
if R < 1, error('dataCell must contain at least one recording.'); end

%-------------------------------
% Check N consistency
%-------------------------------
[N, ~] = size(dataCell{1});
for r = 1:R
    [Nr, ~] = size(dataCell{r});
    if Nr ~= N, error('All recordings must have same N'); end
end

%-------------------------------
% Precompute empirical stats
%-------------------------------
T = zeros(R,1);
dataMean = cell(R,1);
dataCorr = cell(R,1);
empiricalCounts = cell(R,1);
VectorList = mfunc_VectorList(N); % [N x 2^N]
numStates = size(VectorList,2);
for r = 1:R 
    X = dataCell{r}; 
    [~, T(r)] = size(X); 
    dataMean{r} = mean(X,2); 
    dataCorr{r} = (X*X')/T(r); 
    % Empirical counts per state 
    counts_r = zeros(numStates,1); 
    for t = 1:T(r) 
        [~, idx] = ismember(X(:,t).', VectorList.', 'rows'); 
        counts_r(idx) = counts_r(idx) + 1; 
    end
    empiricalCounts{r} = counts_r;
end


T_total = sum(T);
weights = T / T_total; % per-sample normalization

%-------------------------------
% Initialize parameters
%-------------------------------
h = zeros(N,1);
J = zeros(N);
beta = ones(R,1);
h_small_threshold = 1e-4;
J_small_threshold = 1e-3;
%-------------------------------
% Log-likelihood tracking
%-------------------------------
logL_history = zeros(iter_max,1);
prev_logL = -Inf;

%==========================================================================
% Main iteration loop
%==========================================================================
for iter = 1:iter_max
    dh_total = zeros(N,1);
    dJ_total = zeros(N);
    E_emp_vec = zeros(R,1);
    E_model_vec = zeros(R,1);
    
    for r = 1:R
        % Model moments under scaled h,J
        [modelMean_r, modelCorr_r] = mfunc_ModelMeanCorrelation(beta(r)*h, beta(r)*J);
        
        % Gradient contributions
        dh_r = beta(r)*(dataMean{r} - modelMean_r);
        dJ_r = beta(r)*(dataCorr{r} - modelCorr_r);
        dJ_r = dJ_r - diag(diag(dJ_r));
        
        dh_total = dh_total + weights(r)*dh_r;
        dJ_total = dJ_total + weights(r)*dJ_r;
        
        % Energy averages for beta update
        E_emp_vec(r) = -h' * dataMean{r} - 0.5*sum(sum(J.*dataCorr{r}));
        E_model_vec(r) = -h' * modelMean_r - 0.5*sum(sum(J.*modelCorr_r));
    end
    
    % Apply L2 shrinkage to small h
    small_h_idx = abs(h) < h_small_threshold;
    dh_total(small_h_idx) = dh_total(small_h_idx) - 2*lambda_h*h(small_h_idx);
    small_J_idx = abs(J) < J_small_threshold;
    dJ_total(small_J_idx) = dJ_total(small_J_idx) - 2*lambda_J*J(small_J_idx);
       
    % Update h and J
    h = h + lr_rate*dh_total;
    J = J + lr_rate*dJ_total;
    J = (J + J.')/2;
    J(1:N+1:end) = 0; % zero diagonal
    
    % Update beta
    for r = 1:R
        beta(r) = max(0, beta(r) + eta_beta*(E_model_vec(r) - E_emp_vec(r)));
    end
    
    % Compute log-likelihood every check_interval
    if mod(iter, check_interval)==0 || iter==1
        logL_total = 0;
        for r = 1:R
            P_r = mfunc_StateProb(beta(r)*h, beta(r)*J);
            counts_r = empiricalCounts{r};
            validIdx = counts_r > 0;
            logL_total = logL_total + sum(counts_r(validIdx).*log(P_r(validIdx)));
        end
        logL_history(iter) = logL_total;
        rel_improv = abs((logL_total - prev_logL)/max(abs(prev_logL),eps));
        fprintf('Iter %6d | logL = %.8f | ΔlogL = %.3e\n', iter, logL_total, rel_improv);
        if rel_improv < stopping
            fprintf('Converged at iter %d (ΔlogL < %.1e)\n', iter, stopping);
            break;
        end
        prev_logL = logL_total;
    end
end

%-------------------------------
% Package results
%-------------------------------
result.h_star = h;
result.J_star = J;
result.beta = beta;
result.logL_history = logL_history(1:iter);
result.iter = iter;

end

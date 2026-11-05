% 2016/11/11 by T. Ezaki
% This function infers h and J by maximum likelihood estimation.
% binarizedData: nodeNumber x tmax
%
% 2025/9/6 by S. V. Rajesh Modified to accommodate changes
% 2025/10/10: Added automatic L2 shrinkage for small h
%
% DESCRIPTION:
%   Performs iterative maximum likelihood estimation of the Ising model
%   parameters h (node biases) and J (pairwise couplings).
%   Updates h and J using gradient ascent on the likelihood:
%       dh = eta * (dataMean - modelMean)
%       dJ = eta * (dataCorrelation - modelCorrelation)
%
%   NEW FEATURE: Automatic L2 shrinkage of h
%       - For small h (|h_i| < 1e-4), a shrinkage term pulls h_i toward zero
%       - Implemented as: dh_i = dh_i - 2 * lambda_h * h_i
%       - Helps remove meaningless tiny h values arising from near-zero mean data

function result = pfunc_02_Inferrer_ML(binarizedData, lr_rate, iter_max, stopping)

if nargin < 2
    lr_rate = 0.1;
end
if nargin < 3
    iter_max = 1e5;
end
if nargin < 4
    stopping = 5e-6;
end

[nodeNumber,dataLength] = size(binarizedData);

% Empirical statistics
dataMean = mean(binarizedData,2);
dataCorrelation = (binarizedData*binarizedData')/dataLength;

% Initialize parameters
h = zeros(nodeNumber,1);
J = zeros(nodeNumber);

% L2 shrinkage settings for small h
h_small_threshold = 1e-4;  % defines small h
J_small_threshold = 1e-3;  % defines small J
lambda_h = 1e-2;           % strength of L2 shrinkage
lambda_J = 1e-2;           % strength of L2 shrinkage on J
for iter = 1:iter_max
    % Model expectations
    [modelMean, modelCorrelation] = mfunc_ModelMeanCorrelation(h,J);
    
    % Gradient updates
    dh = 1 * (dataMean - modelMean);
    dJ = 1 * (dataCorrelation - modelCorrelation);
    
    % Zero out diagonal of J
    dJ = dJ - diag(diag(dJ));
    
    % Apply L2 shrinkage to small h and J
    small_h_idx = abs(h) < h_small_threshold;
    dh(small_h_idx) = dh(small_h_idx) - 2 * lambda_h * h(small_h_idx);
    % Apply L2 shrinkage to small h and J
    small_J_idx = abs(J) < J_small_threshold;  % Identify small values in J
    dJ(small_J_idx) = dJ(small_J_idx) - 2 * lambda_J * J(small_J_idx);  % Apply shrinkage only to small J values

    % Update parameters
    h = h + lr_rate*dh;
    J = J + lr_rate*dJ;
    
    % Symmetrize J
    J = (J + J.')/2;
    
    % Convergence check
    if (sqrt(norm(dJ,'fro')^2 + norm(dh)^2)/nodeNumber/(nodeNumber+1) < stopping)
        break;
    end
end

fprintf('Converged in %d iterations\n', iter);

% Return results
result.h = h;
result.J = J;

end

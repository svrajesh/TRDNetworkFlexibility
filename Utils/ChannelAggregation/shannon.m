function H = shannon(x)
% SHANNON returns the Shannon entropy for time series x.
% Inputs:
%   x - time series
% Output:
%   H - Shannon entropy of x
    [f, ~] = ksdensity(x); % Estimate probability density function
    f(f == 0) = []; % Remove zeros to avoid log issues
    H = -sum(f .* log2(f)); % Compute Shannon entropy
end

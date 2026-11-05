function bin_data = binarize(data, threshold, sliding)
% BINARIZE Binarizes EEG data using mean thresholding
% Inputs:
%   data    - T x N matrix (T: time, N: channels)
%   winsize    - Window size for sliding binarization
%   sliding - true/false flag for sliding mode
% Output:
%   bin_data - N x T binarized matrix

[T, N] = size(data);  % time x channels

if sliding
    bin_data = zeros(T, N);  % Preallocate logical array

    for ch = 1:N
        % Extract single-channel signal
        x = data(:, ch);

        % Estimate moving baseline (mean)
        baseline = movmean(x, winsize, 'Endpoints', 'shrink');

        % Alternatively, for median baseline use:
        % baseline = medfilt1(x, win_size, 'truncate');

        % Binarize: 1 if signal > baseline, else 0
        bin_data(:, ch) = 2*(x > baseline)-1;
    end

else
    % Static binarization over the full time period
    means = mean(data, 1);
    binthresh = prctile(data,threshold);
    bin_data = 2 * (data > binthresh) - 1;
    bin_data(data == binthresh) = -1;
end

% Return data as channels x time (transpose)
bin_data = bin_data';
end

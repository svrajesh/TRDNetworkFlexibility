function binned_data = bin_downsample(data, fsin, fsout, usePower)
% BIN_DOWNSAMPLE Bins (downsamples) multichannel time series data to lower fs
% Inputs:
%   data    - T x N matrix (time x channels)
%   fsin    - input sampling frequency (Hz)
%   fsout   - output sampling frequency (Hz)
%   usePower - boolean, true: average power (envelope^2), false: average envelope
% Output:
%   binned_data - binned data matrix of size Tbinned x N

if nargin < 4
    usePower = false; % default: envelope
end

[T, N] = size(data);

bin_size = round(fsin / fsout);  % samples per bin

num_bins = floor(T / bin_size);

% Trim data to fit exact bins
data_trim = data(1:(num_bins*bin_size), :);

% Reshape data: bin_size x num_bins x N
data_reshaped = reshape(data_trim, bin_size, num_bins, N);

if usePower
    % Square envelope to get power, then average over bin
    binned_data = squeeze(mean(data_reshaped.^2, 1));
else
    % Average envelope directly
    binned_data = squeeze(mean(data_reshaped, 1));
end

end

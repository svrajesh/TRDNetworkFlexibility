function cleaned = filter_notch(x, FS, lower, upper)
% FILTER_NOTCH notch filters data at specified notch bands using FIR filters
% Applies zero-phase band-stop FIR filtering using filtfilt.
% Inputs:
%   x - time series (1D)
%   FS - data sampling frequency; defaults to 2000 Hz
%   lower - low frequency of each notch band (e.g. [59, 119, 179])
%   upper - high frequency of each notch band (e.g. [61, 121, 181])
% Output:
%   cleaned - notch filtered time series
% Adapted from Author: Alex Gonzalez
% Stanford Memory Lab

if nargin < 2
    FS = 2000;
end

if nargin < 3
    lower = [59, 119, 179];
end

if nargin < 4
    upper = [61, 121, 181];
end

assert(length(lower) == length(upper), 'Lower and upper cutoff arrays must match in length.');

nyquistF = FS / 2;
cleaned = x;

% FIR filter order
N = 500;

%% Loop through each requested notch band
for i = 1:length(lower)
    Wn = [lower(i), upper(i)] / nyquistF; % Normalize frequency band
    b = fir1(N, Wn, 'stop', window(@hann, N+1)); % FIR stopband filter
    cleaned = filtfilt(b, 1, cleaned); % Zero-phase filtering
end

end


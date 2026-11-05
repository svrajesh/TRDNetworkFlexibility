function cleaned = filter_bandpass(x, FS, lower, upper)
% FILTER_BANDPASS bandpass filters data in a specified band using a fixed FIR filter order.
% This function applies a fixed FIR filter order of 500 to bandpass filter the input data.
%
% Inputs:
%   x - time series (signal to be filtered)
%   FS - data sampling frequency; if unspecified, set to 2kHz
%   lower - low frequency of band (lower cutoff)
%   upper - high frequency of band (upper cutoff)
% Output:
%   cleaned - band filtered time series
% Adapted from Author: Alex Gonzalez
% Stanford Memory Lab

if nargin < 2
    FS = 2000; % Default sampling rate of 2000 Hz
end

nyquist = FS / 2;   % Nyquist frequency (half of the sampling rate)

% Normalize the cutoff frequencies (must be between 0 and 1)
Wn = [lower upper] / nyquist;

% Set filter order to a fixed value (500)
N = 500;

% Create the FIR bandpass filter using fir1
b = fir1(N, Wn, 'bandpass', window(@hann, N+1));

% Apply the filter to the signal using filtfilt (zero-phase filtering)
cleaned = filtfilt(b, 1, x);

end



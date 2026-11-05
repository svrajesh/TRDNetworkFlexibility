function envelope = amplitude_envelope(x)
% AMPLITUDE_ENVELOPE returns the envelope of the amplitude of the analytic
% signal using the Hilbert transform. Useful as a measure of spectral
% power. In particular, spectral power is usually the square of the
% magnitude of the Hilbert transform, here we effectively compute the
% square root of the spectral power.

% Inputs:
%   x - time series
% Output:
%   envelope - envelope of data
    s = hilbert(x);
    envelope = abs(s);
end

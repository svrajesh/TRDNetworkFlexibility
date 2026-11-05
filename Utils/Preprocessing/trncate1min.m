function data_out = truncate1min(data, fs, trim_time, keep_time)
% TRUNCATE1MIN - Removes edge artifacts and chops off the first and last
% second (originally used to return the central 1 min of recording but
% adjusted for longer/shorter recs). optionally can supply keep time to
% trim to desired length.
%
% Inputs:
%   data - T x N matrix (timepoints x channels), already filtered & Hilbert-processed
%   fs   - Sampling rate (Hz)
%
% Output:
%   data_out - truncated data (T_trimmed x N)
%
% Assumes data has at least (60 + 1) seconds to allow for edge trimming

        
% Total length check
[T, ~] = size(data);
% Convert to samples
trim_samples = round(trim_time * fs);
if ~isempty(keep_time)
    keep_samples = round(keep_time * fs);
else
    keep_samples = T - 2*trim_samples;
end

usable_start = trim_samples + 1;
usable_end = usable_start + keep_samples - 1;

if usable_end + trim_samples > T
    error('Recording too short after trimming. Need at least %.1f seconds total.', keep_time + 2 * trim_time);
end

% Extract clean 60s segment
data_out = data(usable_start : usable_end, :);
end

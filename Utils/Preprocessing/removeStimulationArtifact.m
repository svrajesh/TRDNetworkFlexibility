function cleanedData = removeStimulationArtifact(data, fs, stimFreq,wide)
% removeStimulationArtifact removes periodic stimulation artifacts from multi-channel data
% by detecting peaks and replacing artifact windows using PCHIP interpolation.
%
%   cleanedData = removeStimulationArtifact(data, fs, stimFreq)
%
%   Inputs:
%       data      - [N x T] matrix (N channels, T samples)
%       fs        - sampling frequency in Hz
%       stimFreq  - stimulation frequency in Hz
%
%   Output:
%       cleanedData - [N x T] matrix with stimulation artifacts interpolated
%
%   Example:
%       cleaned = removeStimulationArtifact(data, 30000, 130);

%% ---------------- Parameters ----------------
preSamples   = 25;     % samples before artifact
if wide
    postSamples  = 100;  
    % samples after artifact
else
    postSamples = 50;
end
contextWin   = 20;     % samples on each side for interpolation context
minPeakHeightFactor = 0.3;
if wide
    minPeakDist  = fs / (2 * max(stimFreq, 1));
else
    minPeakDist = fs/(1.2*max(stimFreq,1))+0;
end

[N, T] = size(data);
cleanedData = data;

%% ---------------- Handle stimFreq == 0 ----------------
if stimFreq == 0
    return
end

%% ---------------- Loop over channels ----------------
for ch = 1:N
    trace = double(data(ch, :));
    
    % --- Detect artifact peaks ---
    % --- Detect artifact peaks (both positive and negative) ---
    
    secder = abs(diff(diff(trace)));
    
    if wide
        threshsecder = double(secder > mean(secder)+std(secder));
        threshsecder = threshsecder(:);                     % ensure column vector
        IBI_min = minPeakDist;               % <- your lower bound on inter-burst interval
        
        % Step 1: find indices where the signal is 1
        idx1 = find(threshsecder == 1);
        
        if isempty(idx1)
            burst_mid = [];
        end
        
        % Step 2: compute gaps between consecutive 1s
        g = diff(idx1);
        
        % A new burst starts whenever the gap exceeds the lower-bound IBI
        burst_start_idx = [1; find(g > IBI_min) + 1];
        burst_end_idx   = [find(g > IBI_min); length(idx1)];
        
        numBursts = length(burst_start_idx);
        
        % Step 3: compute mid-point of each burst using all 1-indices in that burst
        burst_mid = zeros(numBursts,1);
        
        for k = 1:numBursts
            thisBurst = idx1(burst_start_idx(k) : burst_end_idx(k));
            burst_mid(k) = round(mean(thisBurst));
        end
        loc_diff = burst_mid - 20;
        
    else
        threshold = mean(secder)+std(secder); %salvage 
        [pks_diff, loc_diff] = findpeaks(secder,'MinPeakHeight',threshold,'MinPeakDistance', minPeakDist);
    end
    
    % Guard against boundaries
    loc = loc_diff(loc_diff > preSamples & loc_diff < T - postSamples);
    numEvents = numel(loc);
    if numEvents < 5
        continue;
    end
    
    
    % --- Artifact interpolation using PCHIP ---
    for n = 1:numEvents
        range = (loc(n)-preSamples):(loc(n)+postSamples);
        
        % Define context region before and after artifact
        leftCtx  = max(range(1) - contextWin, 1);
        rightCtx = min(range(end) + contextWin, T);
        
        % Known (non-artifact) indices and values
        x_known = [leftCtx:(range(1)-1), (range(end)+1):rightCtx];
        y_known = trace(x_known);
        
        % Interpolate across artifact window
        x_interp = range;
        interpVals = pchip(x_known, y_known, x_interp);
        
        % Replace artifact region with interpolated data
        cleanedData(ch, range) = interpVals;
    end
end



end



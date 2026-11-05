function aperiodics = computeAperiodics(data, fs, method, reg_idx)
%% Compute aperiodic component from data T x K and average by region index
% data: [T x K]
% reg_idx: [1 x K] integers from 1..R
% Output: aperiodics [R x 1]

[~, K] = size(data);
R = max(reg_idx);
aperiodics = zeros(R,1);

switch method
    case 'Hacker'
        win = fs * 1;
        noverlap = round(win * 0.5);
        nfft = max(256, 2^nextpow2(win));
        freqRange = [20 45];
        nLogBins = 50;
        
        slopes = zeros(K,1);
        
        % Compute slope for each channel
        for k = 1:K
            [pxx, f] = pwelch(data(:, k), win, noverlap, nfft, fs);
            idx = f >= freqRange(1) & f <= freqRange(2);
            f_sel = f(idx);
            pxx_sel = pxx(idx);
            
            % Interpolation grid (log-spaced)
            log_f = logspace(log10(freqRange(1)), log10(freqRange(2)), nLogBins);
            
            % Interpolate in log?log space
            log_pxx = 10.^interp1(log10(f_sel), log10(pxx_sel), log10(log_f), 'linear', 'extrap');
            
            % Now take logs for the linear fit
            x = log10(log_f(:));
            y = log10(log_pxx(:));
            
            plot(x, y);
            hold on
            
            p = polyfit(x, y, 1);
            slopes(k) = p(1);
        end
        
        
        % Average by region index
        for r = 1:R
            aperiodics(r) = mean(slopes(reg_idx == r));
        end
end
end

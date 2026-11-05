function [data,analysiswindows] = preprocess(path, wave, stim)
% PREPROCESS2 Preprocesses neural data
% Inputs:
%   path - path to neural data (ns3 or ns5)
%   wave - frequency band to filter in; normal use case is to input 'all'
%   or leave empty, since 'all' is the default, which will then process in all
%   the bands. This can be used standalone to processwithin 1 band if
%   needed in other projects. Doing this way minimizes access to raw nsX
%   data
% Output:
%   data - processed data in struct form

% If band not provided, do 'all'
if nargin < 2
    wave = 'all';
    stim = 0;
elseif nargin < 3
    stim = 0;
end



% Always start with raw loaded data
[rawdata,analysiswindows] = load_data(path,stim);

% Apply notch filter to all channels (skip time col)
for col = 2:width(rawdata)
        rawdata{:, col} = filter_notch([rawdata{:, col}]);
end

%% Loop through all bands if requested
if strcmp(wave, 'all')
    bands = {'raw','full','delta','theta','alpha','beta','gamma','higamma'};
    data = struct;
    for b = 1:numel(bands)
        band = bands{b};
        tmp = rawdata; % copy so original notch-filtered is preserved

        switch band
            case 'delta'
                rng = [0.5 4];
            case 'theta'
                rng = [4 8];
            case 'alpha'
                rng = [8 15];
            case 'beta'
                rng = [15 30];
            case 'gamma'
                rng = [30 50];
            case 'higamma'
                rng = [70 100];
            case 'full'
                rng = [0.5 100];
            case 'raw'
                rng = [0.5 300];
        end

        % Apply bandpass
        for col = 2:width(tmp)
            tmp{:, col} = filter_bandpass([tmp{:, col}], 2000, rng(1), rng(2));
        end

        % Apply envelope unless raw
        if ~strcmp(band, 'raw')
            for col = 2:width(tmp)
                tmp{:, col} = amplitude_envelope([tmp{:, col}]);
            end
        end

        % Store in struct
        data.(band) = tmp;
    end

    %% Execute on the specified band
    % Note to SVR: this could be improved to accept a cell array of bands to
    % process
else
    data = rawdata; % work on notch-filtered base
    switch wave
        case 'delta'
            rng = [0.5 4];
        case 'theta'
            rng = [4 8];
        case 'alpha'
            rng = [8 15];
        case 'beta'
            rng = [15 30];
        case 'gamma'
            rng = [30 50];
        case 'higamma'
            rng = [70 100];
        case 'full'
            rng = [0.5 100];
        case 'raw'
            rng = [0.5 300];
    end

    % Apply bandpass
    for col = 2:width(data)
        data{:, col} = filter_bandpass([data{:, col}], 2000, rng(1), rng(2));
    end

    % Apply envelope unless raw
    if ~strcmp(wave, 'raw')
        for col = 2:width(data)
            data{:, col} = amplitude_envelope([data{:, col}]);
        end
    end
end
end

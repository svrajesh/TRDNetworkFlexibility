function [spectrum,fAll,scores] = speccalc(datafolder)

datadirectory = dir(datafolder);
datadirectory = datadirectory(~startsWith({datadirectory.name},'.'));
spectrum = [];
fAll = [];
scores = zeros(1,length(datadirectory));
parfor df = 1:length(datadirectory)
    meta = get_info(fullfile(datadirectory(df).folder,datadirectory(df).name));
    score = str2num(strrep(extractAfter(meta{3},'cat'),'-','.'));
    scores(df) = score;
    if contains(datadirectory(df).name,'OFF')
        scores(df) = 0;
    elseif contains(datadirectory(df).name,'ON')
        scores(df) = 1;
    end
    % load
    dat = load(betterLoadSave(fullfile(datadirectory(df).folder,datadirectory(df).name)));
    if ~contains(datafolder,'AmpMod')
        data = truncate1min(dat.data, 2000, 1, []); % keep 1 minute of data, get rid of first and last seconds of rec (edge artifacts)
    else
        data = dat.data;
    end

    
    [pxx, f] = pwelch(data, 2000, 1000, [], 2000); % Compute power spectrum using pwelch
    spectrum(df,:,:) = pxx;
    fAll(df,:) = f;
end

end
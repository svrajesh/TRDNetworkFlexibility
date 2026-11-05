function [data,analysiswindows] = load_data(path,stimexpt)
% LOAD_DATA loads neural data. Relies on NPMK Toolbox on path
% Inputs:
%   path - path to neural data (ns3 or ns5)
% Output:
%   data - time series, downsampled to 2kHz if oversampled, stim artifact
%   removed if stimulation present
if nargin == 1
    stimexpt = 0;
end
dat = openNSx(path);
% concatenate stitched data
phantom = 1;
numChunks = length(dat.Data);

if numChunks > 1 && iscell(dat.Data)
    chunkedCat = [];
    for nC = 1:numChunks
        chunkedCat = [chunkedCat, dat.Data{1,nC}];
    end
    dat.Data = chunkedCat;
end
        
elecInfo = dat.ElectrodesInfo;

% include relevant name info for any sEEG channels of interest vv these are
% the PF electrodes of interest
sEEGchans = find(contains({elecInfo.Label},{'PF','OF','AMY','ACC','DLP','AC'},"IgnoreCase",true)); % replace with identifiers for sEEG chans in dataset
colnames = {elecInfo(sEEGchans).Label};

% set the datatable and resampled data to 2kHz if needed
matx = double(dat.Data);
analysiswindows = [];
wide = 0;
if stimexpt
    stimSyncChanInd = find(contains({elecInfo.Label},'StimSync'));
    stimSyncChan = matx(stimSyncChanInd,:);
    %plot(stimSyncChan);
    changepoints = find(abs(diff(stimSyncChan))>1e4);
    analysiswindows = [1 changepoints(1)];
    currend = changepoints(1);
    for i = 1:length(changepoints)-1
        winend1 = changepoints(i);
        winend2 = changepoints(i+1);
        if winend2 - winend1>100e3
            analysiswindows = [analysiswindows; currend+1 winend1; winend1+1 winend2];
            currend = winend2;
        end
    end
    analysiswindows = [analysiswindows; analysiswindows(end,2) changepoints(end); changepoints(end) length(stimSyncChan)];
    bounds = reshape(analysiswindows,1,[]);

    stimvec = [];
    wide = 0;
    % replace stimvec with ordered list of stimulation freqs delivered
    % (ex. [0 130 0 130 0]) and set wide = 0 for when artifacts are
    % visually wide in the time series (ex. SCC stimulation). Adjust
    % analysis windows as appropriate to the experiment. Can vary this for
    % each subject, for example: if contains(pat,subID') stimvec = (stimvec
    % for that subject) wide = (wide for that subject) end, allowing
    % customization across subjects/experiments.


end
matx = matx(sEEGchans,:);

if stimexpt
    for aw = 1:size(analysiswindows,1)
        inds = analysiswindows(aw,1):analysiswindows(aw,2);
        if stimvec(aw) > 0
            matx(:,inds) = removeStimulationArtifact(matx(:,inds), 30000, stimvec(aw), wide);
        elseif phantom
            matx(:,inds) = matx(:,inds);
            % if wanting to test the effect of interpolation on processing, can substitute this function: phantomBlank(matx(:,inds), 30000, 130, wide);
        end
    end
end



if dat.MetaTags.SamplingFreq >2000
    matx = resample(matx',2000,dat.MetaTags.SamplingFreq)';
end
data = array2table(matx');
data.Properties.VariableNames = colnames;

% time axis
tstamps = (1:height(data))/2000; %dat.MetaTags.SamplingFreq;
data.Time = tstamps';
data = [data(:,end) data(:,1:end-1)];

data = bipolar_ref(data);


% Attach participant info as metadata
data.Properties.UserData = get_info(path); %this is largely irrelevant since i save the data as a matrix but could come in handy 

end

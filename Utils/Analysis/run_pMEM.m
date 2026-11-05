function result = run_pMEM(data, threshold, fn,trunc,ablation)
% RUN_PMEM is the second wrapper for fitting PMEM. It passes data
% and window information to a subfunction PFUNC_02_INFERRER_ML, then passes the
% outputs to an accuracy function PFUNC_03_ACCURACY, then stores the data in a struct and
% passes it back to the superfunction

% importantly there is a binning step for smoothing the envelope data so
% that adjacent samples are temporally decorrelated to increase pMEM
% accuracy. Currently setting a 10ms bin.
%
% as of 09082025 we are using a 60s window (global binarization) without
% sliding. Theres only 1 window for the analysis. 
%
% as of 09112025 doing varying binarization thresholds

% Inputs:
%   data - filtered + Hilbert transformed + amplitude envelope data sampled
%   at 2 kHz
%   window - relevant binarization window in MILLISECONDS
%   fn - filename used to build save path for binarized data

% Outputs:
%   result - struct with all the data 

binwidth = .010; % 10 ms
fsout = 1/binwidth;
switch nargin
    case 3
        trunc = 1;
        ablation = 0;
    case 4
        trunc = trunc;
        ablation = 0;
    case 5
        trunc = trunc;
        ablation = ablation;
end

%% Up Front Data Preprocessing
if trunc
    data = truncate1min(data, 2000, 1, []); % keep all of the data, get rid of first and last seconds of rec (edge artifacts)
else
    data = data;
end
data = bin_downsample(data, 2000, fsout); % bin the amp envelope data at fsout
threshUse = round(threshold);%round(winsize*fsout/1000); % winsize to use will be rounded slightly, but thats fine
df = binarize(data, threshold, 0);  % Binarize the data using sliding window or chunked

%% Compute PMEM
params = pfunc_02_Inferrer_ML(df);
[probN, prob1, prob2, rD, r, D1, D2] = pfunc_03_Accuracy(params.h, params.J, df);
mets.d = rD;
mets.KLD1 = D1;
mets.KLD2 = D2;
mets.r = r;
R2imem = corrcoef(probN(:), prob1(:)); 
R2imem = R2imem(1,2)^2;
R2pmem = corrcoef(probN(:), prob2(:)); 
R2pmem = R2pmem(1,2)^2;
mets.R2imem = R2imem; 
mets.R2pmem = R2pmem;
result = struct('filename',fn,'parameters', params, 'metrics', mets);  % Return results in a struct


%% Save the binarized data memory efficiently
binarydata = bitpack_helper(df,'encode');
if ablation
    fnsv = betterLoadSave(strrep(strrep(fn,'Data','BinarizedData/Ablation'),'clean',['Reg_',num2str(ablation),'/binarized']));
else
   fnsv = betterLoadSave(strrep(strrep(fn,'Data','BinarizedData'),'clean',[num2str(threshUse),'/binarized']));
 
end
[fld,~] = fileparts(fnsv);
if ~isfolder(fld)
    mkdir(fld);
end
save(fnsv,"binarydata");

end


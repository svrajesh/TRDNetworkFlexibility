%% Step04 Aperiodic Componnt
% analogous to Hacker first, fit slope between 20 and 45 in log log space
p = gcp("nocreate");
if isempty(p)
    parpool(32);
end
rawdatadir = betterdir(betterLoadSave([USER_DIR,'/Data/',task,'/',pt,'/clean_raw/']));
load(betterLoadSave(['USER_DIR,'/Info/',pt,'_Regions.mat']))
keepreg = [1,3,4,5,7,9,10,11]; %[1,2,4,5,6,8,9,11,12,13];



aperiodicStruct = struct;
for rdf = 1:length(rawdatadir)
    % Load the data for the current raw data file
    rdfn = betterLoadSave(fullfile(rawdatadir(rdf).folder, rawdatadir(rdf).name));
    rawdata = load(rdfn);
    chanNames = rawdata.chanNames;
    chanNames = chanNames(2:end);
    orderingvec = [];
    passchanidxvec = [];
    for kr = 1:length(keepreg)
        
        thisregchans = find(region_labels == keepreg(kr));
        thisregchannames = strrep(region_names(thisregchans),"'","");
        thisregchanidx = find(ismember(chanNames,thisregchannames));
        orderingvec = [orderingvec;thisregchanidx'];
        passchanidxvec = [passchanidxvec,kr*ones(1,length(find(region_labels == keepreg(kr))))];
    end
    data = rawdata.data(:,orderingvec);
    apers = computeAperiodics(data,2000,'Hacker',passchanidxvec);
    aperiodicStruct(rdf).filename = rdfn;
    aperiodicStruct(rdf).aperiodics = apers;
end

%% Fill into params

thresholds = 50;%round([25,30,35,37.5,40,45,50,55,60,62.5,65,70,75]);
regs = 1:8;
bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};
% pmem feature folder and binary data dir
% Loop thresh and band

for thresh = 1:length(thresholds)
    for ablation = 1:length(regs)
        for b = 1:length(bands)
            load(betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/pMEM/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_result.mat']));
            for rec = 1:length(results)
                truekeeps = regs;
                truekeeps(ablation) = [];
                currrecfn = results(rec).filename;
                rawfile = strrep(strrep(currrecfn,bands{b},'raw'),'WeightedReconData','Data');
                recidx = find(strcmp({aperiodicStruct.filename},rawfile));
                results(rec).parameters.aperiodicSlopes = aperiodicStruct(recidx).aperiodics(truekeeps);
            end
            save(betterLoadSave([USER_DIR,'/Analysis/Ablation/',task,'/',pt,'/pMEM/reg_',num2str(ablation),'_pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_result.mat']),'results','-v7.3');
        end
    end
end



p = gcp("nocreate");
if isempty(p)
    parpool(32);
end
rawdatadir = betterdir(betterLoadSave(['USER_DIR/Data/',task,'/',pt,'/clean_raw/']));
load(betterLoadSave(['USER_DIR/Info/',pt,'_Regions.mat']))
keepreg = [1,3,4,5,7,9,10,11];

aperiodicStruct = struct;
for rdf = 1:length(rawdatadir)

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

thresholds = round([25,30,35,37.5,40,45,50,55,60,62.5,65,70,75]);
bands = {'delta','theta','alpha','beta','gamma','higamma','full','raw'};

for thresh = 1:length(thresholds)
    for b = 1:length(bands)
        load(betterLoadSave(['USER_DIR/Analysis/',task,'/',pt,'/pMEM/pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_result.mat']));
        for rec = 1:length(results)
            currrecfn = results(rec).filename;
            rawfile = strrep(strrep(currrecfn,bands{b},'raw'),'WeightedReconData','Data');
            recidx = find(strcmp({aperiodicStruct.filename},rawfile));
            results(rec).parameters.aperiodicSlopes = aperiodicStruct(recidx).aperiodics;
        end
        save(betterLoadSave(['USER_DIR/Analysis/',task,'/',pt,'/pMEM/pMEM_features_',num2str(thresholds(thresh)),'/',bands{b},'_result.mat']),'results','-v7.3');
    end
end

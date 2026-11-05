method = 'weighted';

p = gcp("nocreate");
if isempty(p)
    parpool(32);
end

rootpath = betterLoadSave('USER_DIR/');
switch method
    case 'entropy'
        load(betterLoadSave(['USER_DIR/ChannelSelection/','Resting','/',pt,'/final.mat']));
    case {'supervised','pca','weighted'}
        load(betterLoadSave(['USER_DIR/ChannelSelection/','Resting','/',pt,'/final_recon.mat']));
end

allreg = 1:12;
keepreg = [1,3,4,5,7,9,10,11];
[~,keepreg] = intersect(allreg,keepreg);
for b = 1:length(final_selections)
    final_selections(b).channels = final_selections(b).channels(keepreg,:);
end

thresholds = [25,30,35,37.5,40,45,50,55,60,62.5,65,70,75];
for t = 1:length(thresholds)
    pMEM_model_fits(rootpath,pt,task,thresholds(t),final_selections,method,0);
end

p = gcp("nocreate");
if isempty(p)
    parpool(32);
end

startIndex = 1;
pmemfiles = dir(betterLoadSave(['USER_DIR/Analysis/',task,'/',pt,'/pMEM/*/*.mat']));
pmemfiles = pmemfiles(~startsWith({pmemfiles.name},'.'));
pmemfiles = pmemfiles(~endsWith({pmemfiles.name},'Accuracies.mat'));
Nfiles = length(pmemfiles);
for fl = startIndex:Nfiles

    flname = fullfile(pmemfiles(fl).folder,pmemfiles(fl).name);
    analyze_energy_landscape(flname);

end

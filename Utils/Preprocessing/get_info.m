function info = get_info(filename)
% GET_INFO parses the datafilename and extracts relevant meta data/info
% Inputs:
%   filename - string path to ns3 or ns5 file with the neural data

% Output:
%   info - struct with relevant info including subject name, cat-di score
%   etc

parts = strsplit(filename, '/');
fname = parts{end};
fparts = strsplit(fname,'_'); % the filename inside the folder contains the information

if startsWith(fname,'DBS')
    info = fparts;
else
    % Extract subID, EMU num, CAT-DI score
    try
    ID = extractAfter(fparts{find(contains(fparts,'subj'))},'-'); % this should be the same as pt
    catch
        disp(fparts);
    end
    % CAT-DI scores are stored below in an excel sheet containing a row for each filename of NSX data and a Score column containing the CAT-DI score associated with it
    catfilename = betterLoadSave(sprintf('/USER_DIR/CATDI/%s_CATDI_Scores.xlsx',ID));
    
    try
        catfile = readtable(catfilename);
        if contains(fname,{'CATDI','CAT-DI'})
            idx = find(strcmp(catfile.CATDI,string(fname)));
            task = 'CATDI';
            catdiscore = catfile.Score(idx);
        elseif contains(fname,{'EyesOpen','Baseline'})
            idx = find(strcmp(catfile.Resting,string(fname)));
            task = 'Resting';
            catdiscore = catfile.Score(idx);
        elseif contains(fname,{'Amp','Mod'},"IgnoreCase",true)
            task = 'AmpMod';
            catdiscore = 0;
        end
    catch
        catdiscore = 0; % if the cat score cannot be read assign 0 and deal with it later
        task = '';
        disp(fname);
    end
    epoch = fparts{find(contains(fparts,'EMU'))};
    
    
    % Return information as a cell array
    info = {ID, 'TRD', catdiscore, epoch,task};
end
end

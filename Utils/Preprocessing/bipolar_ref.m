function newdata = bipolar_ref(data)
% BIPOLAR_REF rereferences data in a bipolar scheme on each electrode
% Inputs:
%   data - original set of data
% Output:
%   newdata - rereferenced data

%% Subdivide the table into tables for each electrode
contactnames = data.Properties.VariableNames;
contactnames = strrep(contactnames,'_','-');
timeind = find(strcmp(contactnames,'Time'));
elecnames = unique(regexp(contactnames,'^\D+','match','once')); 
elecnames = elecnames(startsWith(elecnames,'L','IgnoreCase',true) | startsWith(elecnames,'R','IgnoreCase',true)); %isolates each sEEG electrode
newdata = data(:,timeind); % initialize

%% Loop the sEEG elecs
for elec = 1:length(elecnames)
    thiselecinds = find(contains(contactnames,elecnames{elec}));
    thiselecoldtable = data(:,thiselecinds);
    thiseleccontacts = thiselecoldtable.Properties.VariableNames;
    newlabels = {};
    newdat = [];

    % BP reference scheme of n+1 - n
    for cont = 2:length(thiseleccontacts)
        newlabels{cont-1} = strcat(thiseleccontacts{cont},'_',thiseleccontacts{cont-1});
        newdat(:,cont-1) = thiselecoldtable{:,cont} - thiselecoldtable{:,cont-1};
    end
    thiselecnewtable = array2table(newdat);
    thiselecnewtable.Properties.VariableNames = newlabels;

    % concatenate across electrodes
    newdata = horzcat(newdata,thiselecnewtable);
end




function fpn = betterLoadSave(fp)
%% Manages file path differenes between lamella mounted on macOS vs lamella on Linux BioHPC

if ismac
    if contains (fp,'work')
        fpn = strrep(fp,'endosome','Volumes');
    elseif contains(fp,'project') & ~startsWith(fp,'/Volumes')
        fpn = strcat('/Volumes/',fp);
    elseif contains(fp,'archive') & ~startsWith(fp,'/Volumes')
        fpn = strcat('/Volumes/',fp);
    else
        fpn = fp;
    end
elseif isunix & ~ismac
    if contains (fp,'work')
        fpn = strrep(fp,'Volumes','endosome');
    elseif contains(fp,'project')
        fpn = strrep(fp,'/Volumes','');
    elseif contains(fp,'archive') 
        fpn = strrep(fp,'/Volumes','');
    else
        fpn = fp;
    end
elseif ispc & ~(ismac | isunix)
    if contains (fp,'work')
        fpn = strrep(fp,'/Volumes','Y:');
    elseif contains(fp,'project')
        fpn = strrep(fp,'/Volumes','X');
    elseif contains(fp,'archive') 
        fpn = strrep(fp,'/Volumes','Z:');
    else
        fpn = fp;
    end
end

end



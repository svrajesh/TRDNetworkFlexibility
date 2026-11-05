function cleandir = betterdir(foldername)
% function to only extract files, omit '.' files
tempdir = dir(foldername);
cleandir = tempdir(~startsWith({tempdir.name}, '.'));
end
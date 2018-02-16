function [ output_args ] = package(mfile)
%package - gather all dependent files of the main one, and create a single
%folder containing the necessary files

[fList] = chkDepend(mfile)

% New location
locn = ['/Users/jmckenzi/Desktop/' mfile(1:end-2)];
if ~exist(locn,'dir')
    mkdir(locn);
end

for n = 1:numel(fList)
    
    copyfile(fList{n},locn);
end


end


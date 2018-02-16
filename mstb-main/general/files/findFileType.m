function [hits] = findFileType(locn, type)
% findFileType - for a given folder find all files in it with the specified
% extension. This is useful for folders which haven't been indexed, like
% Ottmar's hard drive or the Z drive.

% if strcmp(locn(end),filesep)
%     locn = locn(1:end-1);
% end

% Can we do it purely iteratively and dynamically, rather than going
% through all folders at the beginning...
%[hits] = findFiles(locn,type,{'Path','Name'});
[hits] = fileFinderAll(locn,type);
%clc;

% Sort, trim and display
hits = hits(2:end,:);
%disp(hits);

return

% Now copy these files to a special folder...
loc = '/Users/jmckenzi/DB/OCT/';
for n = 1:size(hits,1)
    
    % Find the name of the folder in which the file is placed
    path = hits{n,1};
    if strcmp(path(end),filesep)
        path = path(1:end-1);
    end
    slash = strfind(path,filesep);
    subfold = path(slash(end)+1:end)
    
    newlocn = [loc subfold filesep hits{n,2}];
    if ~exist([loc subfold],'dir')
        mkdir([loc subfold]);
    end
    
    copyfile([hits{n,1} filesep hits{n,2}],newlocn);
end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hits] = findFiles(fold,type,hits)
% Find files of type in the current folder

% Start by generating the files within the 'root' folder
[root] = genFiles(fold);
numF = size(root,1);

%numH = size(hits,1);

for n = 1:numF
    
    % Check if it is a folder or not
    if root(n).isdir && ~strcmp(root(n).name(1),'.')
        
        % This is the folder, so we should continue to find things in it
        disp([fold filesep root(n).name]);
        [hits] = findFiles([fold filesep root(n).name],type,hits);
        
    elseif length(root(n).name) > length(type)+1
        
        % Then not a folder, so we ask if it is a file of correct type?
        idx = strfind(root(n).name,['.' type]);
        if ~isempty(idx)
            hits = cat(1,hits,{fold,root(n).name});
        end          
            
        
    end
end
        
%hits = [];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all] = genFiles(fold)
% Generate a list of files in 'fold'

all = dir(fold);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function junk

% Get the file list
if ismac
    fL = ['x:' genpath(locn)];
    seps = strfind(fL,':');
else
    fL = ['x;' genpath(locn)];
    seps = strfind(fL,';');
end

% Structure into which to save the file names
files = struct('path',[], 'name',[], 'elem',[],...
    'meta',[], 'exst',[], 'numb', []);
numF = 0;

% Just a cell to store the file names...

% A wait bar for progress...
wb = waitbar(0,'Please wait');

% Parse fL to get the folder structure...
for n = 1:numel(seps)-1
    i1 = seps(n)  +1;
    i2 = seps(n+1)-1;
    path = fL(i1:i2);
    if ~strcmp(path(end),filesep)
        path = [path filesep];
    end
    disp(path)
    
    % Now run through this folder and find the h5 files...
    dls = dir(path);
    
    for r = 1:size(dls,1)
        % Will ignore the folders etc...
        switch dls(r).isdir       
            case 1         
                % Ignore
            otherwise                
                % This is the file extension
                i3 = length(dls(r).name)-(length(type)-1);
                i4 = length(dls(r).name);
                ext = dls(r).name(i3:i4);

                if strcmp(ext,type)
                    try
                        ftest = h5info([path dls(r).name], '/mz');
                        pass = 1;
                    catch err
                        pass = 0;
                    end
                    
                    if pass == 1
                        % Can use this file
                        numF = numF + 1;
                        files(numF).path    = path;
                        files(numF).name    = dls(r).name;
                        fNames{numF,1}      = dls(r).name;
                    end
                end
        end
    end            
end

return

% This is for extracting stuff like metadata and other things

% Now with the files identified, they can be opened and scrubbed...
for n = 1:numF    
    %disp([files(n).path files(n).name]);
    waitbar(n/numF,wb,files(n).name);
    
    % The h5ReadData function doesn't work for attributes and hence
    % metadata - this is a new function to get all of the metadata
    try
        meta = h5info([files(n).path files(n).name], '/meta');
        [files(n).meta] = parseMeta(meta.Attributes,files(n).path,files(n).name);

        % Get the elements in each of the files
        info = h5info([files(n).path files(n).name]);
        elem = info.Datasets;
        ele2 = cell(size(elem,1),1);
        for r = 1:size(elem,1)
            ele2{r} = elem(r).Name;
        end
        files(n).elem = ele2;   
        try
            files(n).numb = files(n).meta.no;
        catch
            try
                files(n).numb = files(n).meta.No;
            catch
                files(n).numb = 0;
            end
        end

        files(n).meta.fPath = files(n).path;
        files(n).meta.fName = files(n).name;
    catch
        disp('invalid h5 file');
    end
end
delete(wb);

% Now you need to find the unique elements of the files - this is to be
% eventually applied to meta too, so make it compatible
[unqM,fail]  = findUniqueMetadata(files, 'meta')
if fail == 0    
    [unq]   = findUniqueElements(files, 'elem');
else
    files = [];
    unq = [];
    unqM = [];
    return
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


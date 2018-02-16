function [ files ] = dupfiles
%dupfiles - find files with duplicated names and display their paths etc
%etc. Also would be good to work out if they contain the same information,
%though this might not be possible easily.  File size may be the clue, as
%may creation / modification date if we can access that easily

% Define / get path in which to look
locn = pwd;

% Determine subfolders and separations
if ismac
    fL = ['x:' genpath(locn)];
    seps = strfind(fL,':');
else
    fL = ['x;' genpath(locn)];
    seps = strfind(fL,';');
end

% Cell array to store all the information going...
numF = 10000;
files = cell(numF,4);
i = 0;

% Loop through each of the subfolders in turn
for n = 1:numel(seps)-1
    i1 = seps(n)  +1;
    i2 = seps(n+1)-1;
    path = [fL(i1:i2) filesep];
    
    % Now get the contents of this directory
    dire = dir(path);
    
    % Run through each possible entry...
    for r = 1:size(dire,1)
        
        % Decide if it is an imzML file...
        if ~dire(r).isdir && ~strcmp(dire(r).name(1),'.') ...
                && strcmp(dire(r).name(end-1:end),'.m')
            mfn = dire(r).name(1:end-2);
            
            i = i + 1;
            files(i,1) = {mfn};
            files(i,2) = {path};
            files(i,3) = {dire(r).bytes};
            files(i,4) = {dire(r).datenum};
            
        end
        
    end            

end

% Trim the list
files = files(1:i,:);

% Sort the list, then we can look at duplicates
files = sortrows(files,1);

% Now let's look through the list and display duplicated files
unqN = unique(files(:,1));
numU = numel(unqN);
numD = 0;
for n = 1:numU
    
    fx = strcmp(files(:,1),unqN{n});
    
    if sum(fx) > 1
        numD = numD + 1;
        
        disp('**************************************');
        disp(unqN{n})
        
        fy = find(fx);
        for r = 1:numel(fy)
            disp(files(fy(r),2));
        end
        
        disp('--------------------------------------');
        disp(char(10));
        
    end
    


end

disp(['Number of duplicate files = ' num2str(numD)]);

% Not conflicted...

end


function [ output_args ] = copyFilesFromToMirror
%copyFilesFromToMirror - copy Olivia's files.

% Not for amateurs
if ~strcmp(getUser,'jmckenzi')
    return
end

% Source location
orig = '/Volumes/Data/Lab Data/Endometrial/Olivia/Olivia DESI Raw data/1 scan per sec/';

% Mirror location
mirror = '/Volumes/JSM/DB/Endometrial/Olivia/';

% Get a list of folders (not subfolders)
folds = folderList(orig);
numF = numel(folds);

% Loop through each one and copy to mirror
for n = 1:numF
    
    % Subfolder
    sub = [orig folds{n} filesep];
    
    % Find mat files
    tmpF = dir(sub);
    numF = size(tmpF,1);
    
    % Loop through, moving mat files...
    for r = 1:numF
        
        if ~tmpF(r).isdir
            
            % Check for .mat extension
            flag = strcmp(tmpF(r).name(end-2:end),'mat');
            
            % Is it POS or NEG?
            if flag
                isPos = any(strfind(lower(tmpF(r).name),'pos'));
                isNeg = any(strfind(lower(tmpF(r).name),'neg'));
                
                if isPos & ~isNeg                    
                    copyfile([sub tmpF(r).name],[mirror 'POS/' tmpF(r).name]);
                elseif isNeg & ~isPos
                    copyfile([sub tmpF(r).name],[mirror 'NEG/' tmpF(r).name]);
                end
                
            end
            
        end
    end
    disp(int2str(n));
end

end


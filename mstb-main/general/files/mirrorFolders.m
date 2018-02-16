function [ output_args ] = mirrorFolders(orig,new)
%mirrorFolders - create subfolders from one directory in another

ols = dir(orig)
numF = size(ols,1)

% Loop through
for n = 1:numF
    
    name = ols(n).name;
    if length(name) < 3
        continue;
    end
    
    % Check that end of name is not '.raw'
    lst = name(end-2:end);
    
    if ols(n).isdir && ~strcmp(lst,'raw')
        disp(name);
        
        % Duplicate folder name in new location
        if ~exist([new name],'dir')
            mkdir([new name]);
        end
        
    end

end

end


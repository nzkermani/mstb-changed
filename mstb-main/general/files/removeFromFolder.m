function removeFromFolder
%removeFromFolder - for all subdirectories of the specified directory, 
% move their contents to the directory.

% Directory
direct = '/Users/jmckenzi/DB/LNDatabaseML/';

all = dir(direct);
numE = size(all,1);

% Loop through
for n = 1:numE
    
    % Compose
    tmp = [direct all(n).name];
    
    if isdir(tmp) && length(all(n).name) > 3;
        
        % Get contents
        cont = dir(tmp);
        numF = size(cont,1);
        
        % Move contents
        for r = 1:numF
            
            name = cont(r).name;
            
            if length(name) > 3
                movefile([tmp filesep name],[direct name]);
            end
            
        end
        
    end
    
end


end


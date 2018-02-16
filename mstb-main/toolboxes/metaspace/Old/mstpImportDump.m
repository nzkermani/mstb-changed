function [ info ] = mstpImportDump
%mtspImportDump - read in the TXT files dumped from the engine

tic

% What is the path of the files
path = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/Annos/';

% Find all the text files
%allF = fileFinderAll(path,'txt',true);

dd = dir(path);

%allN = {dd.name}';
numF = size(dd,1);

% File info...
info = cell(numF,5);

% Skip the first 2
for n = 1:numF
    
    %try
        % If all else fails...
        
        % Find underscores
        und = strfind(dd(n).name,'_');
        if numel(und) < 3
            continue;
        end
        
        % If there are more underscores than 3, then these are special LNTO
        % files which need to be treated slightly differently
        if numel(und) > 3
            info{n,2} = dd(n).name(1:und(3)-1);
            info{n,3} = dd(n).name(und(3)+1:und(4)-1);
            info{n,4} = dd(n).name(und(4)+1:und(5)-1);
            info{n,5} = dd(n).name(und(5)+1:end-4);           
            
        elseif numel(und) == 3
            info{n,2} = dd(n).name(1:und(1)-1);
            info{n,3} = dd(n).name(und(1)+1:und(2)-1);
            info{n,4} = dd(n).name(und(2)+1:und(3)-1);
            info{n,5} = dd(n).name(und(3)+1:end-4);
            
        else
            % Unknown number of und
            disp('FAIL>>>');
        end
                
        % Split up...
        info{n,1} = dd(n).name;

    %catch
        
        disp(dd(n).name);
                
        
    %end
    
end
   
toc

% Remove first two rows
%info = info(3:end,:);
    

end


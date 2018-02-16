function [info] = jlaSampleMatch
%jlaSampleMatch - DESI file matching for reprocessing James's DESI samples

% Main h5 folder
h5Main = ['/Volumes/Data/Lab Data/Colorectal/JLA Samples/'...
    'Fresh frozen/TO BE ORGANISED/Exactive/'...
    'DESI -ve mode update 17.7.17/H5 files/'];
h5M = fileFinderAll(h5Main,'h5',true);

% Secondary place...
h5Aux = ['/Volumes/Data/Lab Data/Colorectal/JLA Samples/'...
    'Fresh frozen/TO BE ORGANISED/Exactive/'...
    'DESI -ve mode update 17.7.17/Samples/'];
h5A = fileFinderAll(h5Aux,'h5',true);

% Get a list of folder names
fold = folderList(h5Aux);
numF = size(fold,1)
info = cell(numF,4);

% Now we just need to match the folder names (which contain an imzML file)
% to the h5 file names...
for n = 1:numF
    
    % Which folder is this?
    info{n,1} = [h5Aux fold{n}];
    %disp(fold{n});
    
    % Find the imzML file
    tmpI = fileFinderAll(info{n,1},'imzML',true);
    if size(tmpI,1) == 1
        info{n,2} = tmpI{1,2};
    else
        disp(fold{n});
        disp('Too many imzML');
        continue;
    end
    
    % Is there a suitable h5 file in h5Main (h5M)?
    tmpM = ~cellfun(@isempty,strfind(h5M(:,2),fold{n}));
    
    % If so, then let's use that...
    if sum(tmpM) == 1        
        % Save details...
        info{n,3} = h5M{tmpM,1};
        info{n,4} = h5M{tmpM,2};
    
    elseif sum(tmpM) > 1
        % Then there are more than one, so we need to decide which to use.
        % Probably go for the one with the longer name, as it is a modified
        % version of the original one and so likely to be better
        ll = cellfun(@length,h5M(tmpM,2));
        [~,b] = max(ll);
        tmpM = find(tmpM);
        info{n,3} = h5M{tmpM(b),1};
        info{n,4} = h5M{tmpM(b),2};
        
    elseif sum(tmpM) == 0
        % Then we need to look in a different place for an H5 file. This is
        % likely to be the same folder...
        
        newF = fileFinderAll(info{n,1},'h5',true);
        if size(newF,1) == 1
            % Success, just take it...
            info{n,3} = newF{1,1};
            info{n,4} = newF{1,2};
                        
        elseif size(newF,1) > 1
            ll = cellfun(@length,newF(:,2));
            [~,b] = max(ll);
            info{n,3} = h5M{b,1};
            info{n,4} = h5M{b,2};        
        
        else
            disp(fold{n});
            disp('No H5 file found');

        end
        
    else
        disp(fold{n});
        disp('Something wrong');
    end
    
    %info(1:n,:)
    
    
end

end


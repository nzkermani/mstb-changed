function mirrorRawFiles(zfold,dest)
%mirrorRawFiles - copy files from the Z drive (certain types of raw file)
%into the approp folder on the mirrored folder

mode = 'pos';

% Find all raw folders in the Z drive
allZ = fileFinderAll(zfold,'raw',true);

% What about files that we already have locally?
locRaw = 'E:\Data\Olivia\Raw\';
rawF = fileFinderAll(locRaw,'raw',true);

% Now we need to strip out the ones we don't want
switch mode
    case 'neg'       
        f1 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'NEG'));
        f2 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'neg'));
        f3 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'CORR'));
        f4 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'corr'));
        
    case 'pos'
        f1 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'POS'));
        f2 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'pos'));
        f3 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'CORR'));
        f4 = ~cellfun(@isempty,strfind(lower(allZ(:,2)),'corr'));        
end

% Combine
fx = (f1 | f2);
fy = (f3 | f4);
fz = fx & fy;
allZ = allZ(fz,:);
numZ = size(allZ,1);

% Loop through the files
for n = 1:numZ
    
    % Determine the new file location...
    sl = strfind(allZ{n,1},filesep);
    subFold = allZ{n,1}(sl(end)+1:end);
    newLocn = [dest subFold filesep allZ{n,2}];
    
    if exist(newLocn,'dir')
        continue;
    end
    
    disp(allZ{n,2});
    
    % Can we find the exact file name already in the local folder?
    fx = strcmp(rawF(:,2),allZ{n,2});
    if sum(fx) == 1
        % Then we have a direct match, so just move the raw file
        movefile([rawF{fx,1} filesep rawF{fx,2}],newLocn);
        
        disp('---> Moved from local');
        
    else
        
        % Copy from the Z drive...
        copyfile([allZ{n,1} filesep allZ{n,2}],newLocn);
        disp('---> Z --> E');
        
    end
        
        
    
end




end


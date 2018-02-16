function [ output_args ] = mtspImportFiles(info)
%mtspImportFiles - use the cell array of files and slowly read them in...

% Root path?
root = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/Annos/';
savePath = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/MTSP/';

% Remove empty file names
empt = cellfun(@isempty,info(:,1));
info = info(~empt,:);
%numF = size(info,1);

% Find unique images
[unqI,~,indI] = unique(info(:,2));
numI = numel(unqI);

% Range of potential sizes...
allSz = 30:1:150;

wb = waitbar(0,'Importing');

% Loop through each file...
for n = 1:numI
    
    % Generate file name...
    saveName = [savePath unqI{n} '.mat'];
    disp([int2str(n) '/' int2str(numI) ' - ' saveName]);
    if exist(saveName,'file')
        continue;
    end
    
    % Find indices of txt files for this image
    fx = find(indI == n);
    numT = numel(fx);
    
    % Change the waitbar
    waitbar(0,wb,unqI{n});
    
    % Start looping through the txt files
    for r = 1:numT
        
        % Read in
        tmp = dlmread([root info{fx(r),1}],',');
        
        % Determine its index, which is the 2nd parameter
        indx = info{fx(r),3};
        indx = str2double(indx) + 1;
        
        % If first file, then create empty data
        if r == 1
            mz = zeros(numT,1);
            sp = zeros(numel(tmp),numT);
            anno = cell(numT,2);
        end
        
        % Add to the matrix
        mz(indx,1) = indx;
        sp(:,indx) = tmp';
        anno{indx,1} = info{fx(r),4};
        anno{indx,2} = info{fx(r),5};
        
        waitbar(r/numT,wb);
        
    end
    
    % Can we try to determine the best size for reshaping?
    sz = size(sp,1);    
    potSz = sz ./ allSz;
    flSz = floor(potSz);
    dfSz = potSz - flSz;
    idxSz = dfSz == 0;
    sz = [allSz(idxSz); sz ./ allSz(idxSz)]';
    
    % Save the file
    save(saveName,'mz','sp','anno','sz')
    
end

% Delete the waitbar
delete(wb);


end


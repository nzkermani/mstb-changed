function [ mz,sp,anno,sz ] = mtspImportConvert(info)
%mtspImportFiles - use the cell array of files and slowly read them in...

% Save path?
%savePath = '/Volumes/JSM/DB/Metaspace/NewEngineDump/Convert/';

% How many files?
numF = size(info,1);

% Range of potential sizes...
allSz = 30:1:150;

wb = waitbar(0,'Importing');

% Loop through each file...
for n = 1:numF
    
    % Generate file name...
    %saveName = [savePath info{n,1} '.mat'];
    %disp([int2str(n) '/' int2str(numF) ' - ' saveName]);
    %if exist(saveName,'file')
        %continue;
    %end

    % Temporarily dump the info
    tFiles = info{n,2};
    tIndex = info{n,3};
    tForms = info{n,4};
    tAdcts = info{n,5};
    anno = [tForms tAdcts];
    
    % How many images does this file have?
    numI = size(tFiles,1);
    
    % Change the waitbar
    waitbar(0,wb,info{n,1});
        
    % Start looping through the txt files
    for r = 1:numI
        
        % Read in
        tmp = dlmread([tFiles{r,1} filesep tFiles{r,2}],',');
        
        % If first file, then create empty data
        if r == 1
            sp = zeros(numel(tmp),numI);            
        end
        
        % Add to the matrix
        sp(:,r) = tmp';
        
        
        waitbar(r/numI,wb);
        
    end
    
    % Just determine the mz
    [mz] = mtspForm2MZ(anno);

    % Can we try to determine the best size for reshaping?
    sz = size(sp,1);    
    potSz = sz ./ allSz;
    flSz = floor(potSz);
    dfSz = potSz - flSz;
    idxSz = dfSz == 0;
    sz = [allSz(idxSz); sz ./ allSz(idxSz)]';
    
    % Save the file
    %save(saveName,'mz','sp','anno','sz')
    
end

% Delete the waitbar
delete(wb);


end


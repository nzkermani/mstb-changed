function mtspWorkflow
%% mtspWorkflow
% Uses the raw text files output from the Metaspace engine via the Python
% code, this function demonstrates the workflow for getting all of the
% files into a useful format...

% What is the tissue type?
type = 'ovarian';

% Now we just need to add the various parts / locations for each tissue
% type before running the function. Descriptions below:
% newDump   - folder with individual txt files from Metaspace
% oldMerge  - folder containing the OLD merged Metaspace files
% saveConv  - folder to save the Metaspace converted files
% matFiles  - where the toolbox mat files with annotations are stored
% imzFiles  - where the recalibrated imzML files are stored
% saveFinal - where to save the final finished files
switch type
    
    case 'breast'
        
        newDump = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Raw/Breast/';        
        oldMerge = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP/Breast/';        
        %saveConv = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Convert/Breast/';        
        matFiles = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP-Merge/Neg/Breast/';        
        imzFiles = '/Volumes/JSM/DB/Breast DESI/imzML Recal/';
        calFile = '/Volumes/JSM/DB/Metaspace/Recalibrations/Breast-Orig-Precal.mat';
        saveFinal = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/Breast/';
        
    case 'lymph'
        
        newDump = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Raw/Lymph/';        
        oldMerge = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP/Lymph/';        
        %saveConv = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Convert/Lymph/';        
        matFiles = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP-Merge/Neg/Lymph/';        
        imzFiles = '/Volumes/JSM/DB/Lymph DESI Recalibrated/';
        calFile = '/Volumes/JSM/DB/Metaspace/Recalibrations/Lymph-Orig-Precal.mat';
        saveFinal = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/Lymph/';
        
    case 'colorectal'
        
        % Just need to do a limited amount of data processing. Typically,
        % just the imzML reprocessing to replace the old DESI variables.
        % Could actually just fake pretent that the metaspace variables
        % have changed by pointing to the same location.  THat way at least
        % the files will have the same structure...
        newDump = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Raw/Colorectal/';        
        oldMerge = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP/Colorectal/';        
        %saveConv = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Convert/Lymph/';        
        matFiles = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP-Merge/Neg/Colorectal/';        
        imzFiles = '/Volumes/JSM/DB/Colorectal DESI/Metaspace/imzML/';
        calFile = '/Volumes/JSM/DB/Metaspace/Recalibrations/Colorectal-Orig-Precal.mat';
        saveFinal = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/Colorectal/';

        % Create folders for the function to work properly
        ffo = fileFinderAll(oldMerge,'mat',true);
        for n = 1:size(ffo,1)
            ttf = [newDump ffo{n,2}(6:end-4)];
            if ~exist(ttf,'dir')
                mkdir(ttf);
            end            
        end
            
    case 'ovarian'
        
        % Just need to do a limited amount of data processing. Typically,
        % just the imzML reprocessing to replace the old DESI variables.
        % Could actually just fake pretent that the metaspace variables
        % have changed by pointing to the same location.  THat way at least
        % the files will have the same structure...
        newDump = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Raw/Ovarian/';        
        oldMerge = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP/Ovarian/';        
        %saveConv = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Convert/Lymph/';        
        matFiles = '/Volumes/JSM/DB/Metaspace/EngineDumpOLD/MTSP-Merge/Neg/Ovarian/';        
        imzFiles = '/Volumes/JSM/DB/Ovarian-Luisa/imzML/';
        calFile = '/Volumes/JSM/DB/Metaspace/Recalibrations/Ovarian-Orig-Precal.mat';
        saveFinal = '/Volumes/JSM/DB/Metaspace/EngineDumpNEW/Final/Ovarian/';

        % Create folders for the function to work properly
        ffo = fileFinderAll(oldMerge,'mat',true);
        for n = 1:size(ffo,1)
            ttf = [newDump ffo{n,2}(1:end-4)];
            if ~exist(ttf,'dir')
                mkdir(ttf);
            end            
        end
        
        
    otherwise
        disp('Unknown tissue type');
        return;
end

% Read in the precalibration file, which contains the ppm deviations of the
% key ions... We want to include this information in the structure
precal = open(calFile);

% DESI processing options
ooo.ext = 'imzML';
[opts,~] = desiGetProcOptions(ooo,true);
desiOpts = opts.imzml;

% Somewhere to store all of the information
info = struct('file',[],...
    'qtyOld',[],...
    'qtyNew',[],...
    'msiOld',[],...
    'msiNew',[],...
    'precal',[],...
    'status',[]);

% Find all of the files that came from metaspace, and then we can begin the
% processing task...
allF = folderList(newDump);
numF = size(allF,1);
for n = 1:numF
    
    % Format the final file's save name
    ffsn = [saveFinal allF{n,1} '.mat'];
    disp(ffsn);
    
    % Check to see if it exists...
    if exist(ffsn,'file')
        disp('---Skip');
        continue;
    end
    
    % Find the corresponding calibration values for this file...
    switch type
        case {'colorectal','ovarian'}
            pcf = [allF{n,1} '.imzML'];
        otherwise
            pcf = [allF{n,1}(1:end-6) '.imzML'];
    end
    idx = strcmp(precal.files(:,1),pcf);
    if sum(idx) == 1
        info(n).precal = precal.devs(:,idx)';
    end
    
    % Open the toolbox-compatible mat file...
    switch type
        case 'colorectal'
            ofs = 9;
            ofs2 = 0;
        case 'ovarian'
            ofs = 0;
            ofs2 = 0;
        otherwise
            ofs = 6;
            ofs2 = ofs;
    end
    matName = [matFiles allF{n,1}(1:end-ofs) '.mat'];
    if exist(matName,'file')
        tmp = open(matName);
        dpn = tmp.dpn;
        dpn.mtspOLD = dpn.mtsp;
        dpn.mtsp = [];
    else
        info(n).status = 'NO TOOLBOX FILE';
        info(n).file = allF{n,1};
        disp('---No toolbox file');
        continue;
    end
    
    % Open the old data set as downloaded from Metaspace - we just want to
    % know how many annotations there were from the first go
    switch type
        case 'ovarian'
            icl = '';
        otherwise
            icl = 'ICL--';
    end
    oldName = [oldMerge icl allF{n,1}(1:end-ofs2) '.mat'];
    if exist(oldName,'file')
        tmp = open(oldName);
        info(n).qtyOld = numel(tmp.mz);
    end

    % Get a list of the annotations
    switch type 
        case 'colorectal'
            
            mts = open([oldMerge 'ICL--' allF{n,1} '.mat']);
            annoList = mts.anno;
            
        case 'ovarian'
            mts = open([oldMerge allF{n,1} '.mat']);
            annoList = mts.anno;
            
            
        otherwise
            annoList = fileFinderAll([newDump allF{n,1}],'txt',true);
            [idx,anno,adct] = extractInfoDump(annoList(:,2));
            
            % Convert into useable data
            [mts.mz,mts.sp,mts.anno,mts.sz] = importInfoDump(annoList,idx,anno,adct);
    end
    
    % Now is time to sort the images by mz rather than another factor
    [~,idx] = sort(mts.mz);
    mts.mz = mts.mz(idx);
    mts.sp = mts.sp(:,idx);
    mts.anno = mts.anno(idx,:);    
    
    % Now we need to put the mtsp annotations into the correct format
    [newDPN] = mtspMatch(dpn,mts);

    % Save the good bit
    dpn.mtsp = newDPN.d1;
    
    % Now we need to find the imzML file and do the reprocessing for the
    % sake of completeness
    imzStruct.dir = imzFiles;
    imzStruct.nam = [allF{n,1} '.imzML'];
    
    [newMZ,newSP,timestamp] = desiFunctionIMZML(imzStruct,desiOpts);
    dpn.date = timestamp;
    
    % Replace the old data and save the new data
    dpn.d1OLD = dpn.d1;
    dpn.d1.mz = newMZ;
    dpn.d1.sp = newSP;
        
    % Save the completed file in the determined location
    save(ffsn,'dpn');
            
    % Save diagnostic information
    info(n).file    = allF{n,1};
    info(n).qtyNew  = size(annoList,1);
    info(n).msiOld  = numel(dpn.d1OLD.mz);
    info(n).msiNew  = numel(dpn.d1.mz);
    info(n).status  = 'Completed';
    
end

assignin('base','info',info);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx,anno,adct] = extractInfoDump(names)
% Extract info from the file names

und = strfind(names,'_');

sz = size(names,1);

idx = zeros(sz,1);
anno = cell(sz,1);
adct = cell(sz,1);

for n = 1:sz
    
    ii = und{n,1}(end-2:end);
    
    idx(n,1) = str2double(names{n,1}(ii(1)+1:ii(2)-1));
    anno{n,1} = names{n,1}(ii(2)+1:ii(3)-1);
    adct{n,1} = names{n,1}(ii(3)+1:end-4);
    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp,anno,sz] = importInfoDump(list,idx,anno,adct)

wb = waitbar(0,'Importing text dump');

% Range of potential sizes...
allSz = 30:1:150;

% How many files are there?
numI = size(list,1);

% Start looping through the txt files
for r = 1:numI

    % Read in
    tmp = dlmread([list{r,1} filesep list{r,2}],',');

    % If first file, then create empty data
    if r == 1
        sp = zeros(numel(tmp),numI);            
    end

    % Add to the matrix
    sp(:,r) = tmp';
    
    waitbar(r/numI,wb);

end

% Just determine the mz
anno = [anno adct];
[mz] = mtspForm2MZ(anno);

% Can we try to determine the best size for reshaping?
sz = size(sp,1);    
potSz = sz ./ allSz;
flSz = floor(potSz);
dfSz = potSz - flSz;
idxSz = dfSz == 0;
sz = [allSz(idxSz); sz ./ allSz(idxSz)]';

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



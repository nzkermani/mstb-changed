function xxxBatchProcess
%xxxBatchProcess - process a lot of files in a way that won't fail!

% List the directories
% dirs.main = 'D:\Raman Brain\';
% dirs.subs = {'2016-09-15 Raman human brain batch 2';...%'2016-09-13 Raman Mouse brain new batch';...
%     '2016-04 Raman Samples'};
% dirs.output = 'D:\Raman Brain\Reproc\';

% dirs.main = 'Z:\Data\Jocelyn Waters Data\Raw data for Raman paper\';
% dirs.subs = {'Human brain'};
% dirs.output = 'E:\Data\DESI Raman\';

dirs.main = 'Z:\Data\Haixing\Laser Imaging\';
dirs.subs = {'Normal';'Tumour'};
dirs.output = 'E:\Data\Haixing\';

% How many subfolders are there?
numS = numel(dirs.subs);

% Get the processing options - same for all files
tmp.ext = 'raw';
[mainOpts,~] = desiGetProcOptions(tmp);

% For each subfolder, run the batch...
for n = 1:numS
    
    try
        batchFolder(dirs.main,dirs.subs{n},dirs.output,mainOpts.raw);
    catch err
        err
        disp('GENERIC FOLDER ERROR');
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = batchFolder(dirMain,dirSub,dirOp,opts)
% Get the file names and place to output the files to

% Find all files in the folder...
fullDir = [dirMain dirSub];
disp(fullDir);
list = fileFinderAll(fullDir,'dat');
list(1,:) = [];

% Just find the lock mass corrected files
%lmf = strfind(list(:,1),'JSM-LMC');
%lmf = ~cellfun(@isempty,lmf);
%list = list(lmf,:);

numF = size(list,1);

% What is the save directory?
saveDir = [dirOp dirSub filesep];
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

for n = 1:numF

    tic
    
    % Format the file names
    file.dir = [list{n,1} filesep];
    file.nam = list{n,2};
    file.ext = fileExtension(file.nam);

    sl = strfind(file.dir,filesep);
    sl = sl(end-1);
    file.nam = file.dir(sl+1:end-1);
    file.dir = file.dir(1:sl);
    file.ext = 'raw';
    
    % Check file existence
    fex = exist([file.dir file.nam],'file');    
    if ~fex
        disp(['---FAIL ' file.nam]);
        continue;
    end
    
    % Message of successs
    disp(['+++OPEN ' file.nam]);
        
    % Prepare the new file name and folder location
    fileSave = [saveDir file.nam(1:end-4)];
    if exist([fileSave '.mat'],'file')
        disp(['???EXIST ' file.nam]);
        continue;
    end
    disp(['>>>SAVE ' fileSave]);
    
    flag = true;
    try
        [mz,x,~,~,numPoints,~] = h5waters([file.dir file.nam],opts);
        %mz = 100:1:110;
        %x = rand(10,10,numel(mz));
        %numPoints = rand(10,10);
        
    catch err
        disp('!!!FAIL - file processing failure');
        err
        flag = false;        
    end
    
    % If processing failed then just skip through
    if ~flag
        continue;
    end
    
    % Create the guidata structure and draw a few things here
    dpn.file = file;
    dpn.opts = opts;

    dpn.d1.mz = mz;
    dpn.d1.sp = x;
    dpn.d1.numPoints = numPoints;

    dpn.fig = [];
    dpn.defP = dpn.file.dir;

    dpn.mode = 'single';

    % Now save the file to where it should be...
    save(fileSave,'dpn','-v7.3');
    toc
    disp(['+++PASS ' file.nam]);
    
    clear dpn x mz;

end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

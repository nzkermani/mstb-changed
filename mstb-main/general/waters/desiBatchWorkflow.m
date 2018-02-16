function [ status,tt ] = desiBatchWorkflow
%desiBatchWorkflow - process a folder's worth of raw files using James' new
%method, rather than that developed in the toolbox

% What is the resolving power?
rp = 30000;

% Folder locations...
rawFold = 'E:\Data\Olivia\Subsets\Pos\';
savFold = 'E:\Data\Olivia\Subsets\Pos-Proc\';

% Get the list of raw files
rawF = fileFinderAll(rawFold,'raw',true);
numF = size(rawF,1);

% Places to save
tt = zeros(5,numF);
status = cell(numF,1);

% Loop through
for n = 1:numF
    
    % File name
    fn = [rawF{n,1} filesep rawF{n,2}];
    disp([int2str(n) ' - ' fn]);        
    %savName = [savFold subFolder(rawF{n,1}) filesep rawF{n,2}(1:end-4) '-NewMethod.mat'];
    savName = [savFold rawF{n,2}(1:end-4) '-NewMethod.mat'];
    
    if exist(savName)
        continue;
    end

    % Read in the file
    [sp,~,~,xy2D] = desiReadRaw(fn,true);
    
    % Run the function...
    [pp,tt(:,n)] = desiReadRawWorkflow(sp,xy2D,rp);
    
    % How to save pp?
    mz = pp.mz;
    sp = pp.data;
    %save(savName,'mz','sp');
    
    % How about save in a dpn-style structure ready to be imported into the
    % toolbox
    
    dpn.file.dir = [rawF{n,1} filesep];
    dpn.file.nam = rawF{n,2};
    
    dpn.d1.mz = mz;
    dpn.d1.sp = sp;
    
    dpn.fig = [];
    dpn.defP = [rawF{n,1} filesep];
    dpn.date = [];    
    dpn.mode = 'single';
    
    save(savName,'dpn');
    
    status{n,1} = 'Pass';
    
end

end


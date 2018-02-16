function jlaSampleReprocess(list)
%jlaSampleReprocess - load the h5 file, convert to new toolbox format,
%reprocess the imzML file, and save as a new set of files ready for
%analysis.

% Save dir...
sd = '/Volumes/JSM/DB/Colorectal DESI/JLA Reprocess/';

% Define the DESI options
[opts,~] = desiGetProcOptions('imzML',true);
opts = opts.imzml;
javaclasspath('packages/imzMLConverter/imzMLConverter.jar');

% How many files (we skip the ones that don't have an h5 file)
numF = size(list,1);
for n = 1:numF
    
    % Check that there is an h5 file (col3) ... (or mat file - col2)
    if isempty(list{n,3})
        disp(['No H5/imzML file - ' list{n,2}]);
        continue;
    end
    
    % Determine the name which we will save the finished file as. If it
    % already exists then we can skip
    saveName = [sd list{n,4}(1:end-2) 'mat'];
    if exist(saveName,'file')
        
        % Can we determine the imzML acquisition date
        getDate = imzmlDate([list{n,1} filesep list{n,2}])
        
        tmp = open(saveName);
        
        % Replace with the actual imzML file information
        tmp.dpn.file.imzML.dir = [list{n,1} filesep];
        tmp.dpn.file.imzML.dir = list{n,2};
        tmp.dpn.date = getDate;
        dpn = tmp.dpn;
        save(saveName,'dpn');
        
        disp(['Exist - ' list{n,4}]);
        continue;
    end
    
    % Import...
    tic
    h5f.dir = [list{n,3} filesep];
    h5f.nam = list{n,4};
    [dpn] = importFromH5(h5f);
    toc
    
    % Now do the imzML processing...
    tic
    imf.dir = [list{n,1} filesep];
    imf.nam = list{n,2};
    [mz,sp] = importFromIMZML(imf,opts);
    toc
    
    clear h5f imf
    
    % Check file size, and then can replace if matching...
    szO = size(dpn.d1.sp);
    szR = size(sp);
    
    if szO(1) == szR(1) && szO(2) == szR(2)
        dpn.d1.mz = mz;
        dpn.d1.sp = sp;
        
        % Save...
        save(saveName,'dpn');
    else
        disp(['Mismatch - ' list{n,4} ' -- ' list{n,2}]);
    end
    
    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dpn] = importFromH5(file)
% Import the data from the H5 file

% This function actually does the bulk of the work, and can be called
% independently of the GUI
[mz,x,img,opt,allAnno,imzML] = desiH5LoadFunction(file);

% Create a guidata structure, and then save the results, and then display
% them in the right boxes...
dpn.file = file;
dpn.file.imzML = imzML;

% This bit is a little unknown at the moment. Presumably they need to be
% specific to this toolbox...
dpn.opts = [];

% Copy the mz and sp vector/matrices
dpn.d1.mz = mz;
dpn.d1.sp = x;
dpn.d1.img = img;

% Coregistered images...
dpn.opt.coreg = opt.coreg;

% These were provided
%dpn.fig = fig;
dpn.defP = file.dir;

% And the annotations
dpn.anno = allAnno;

% Mandatory for single non PS data
dpn.mode = 'single';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp] = importFromIMZML(file,opts)
% Import imzML file
    
% Run the main function, to get the matrix and mz vector
[mz,sp] = desiFunctionIMZML(file,opts);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

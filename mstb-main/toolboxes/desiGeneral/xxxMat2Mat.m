function xxxMat2Mat
%xxxMat2Mat - replace a processed file (with annotations) with the data
%from another. And do this for a whole folder's worth of files.
%
% TESTED
%
% James McKenzie, 2017.

% Ask the user for a folder, but if nothing found then quit
path = 'E:\Data\Haixing\Haixing-Annotated\Tumour\';
%[path] = uigetdir(pwd);
if path == 0
    return
end

% Determine the output path - just append something in the middle...
opPath = [path 'Mat2Mat\']
if ~exist(opPath,'dir')
    mkdir(opPath);
end

% Find all of the files
list = fileFinderAll(path,'mat');
list(1,:) = [];
numF = size(list,1);

% Need a variable to include the name of the files and what happened to
% them
status = cell(numF,4);

% Get the folder from the user
%[reProcPath] = uigetdir(pwd);
reProcPath = 'E:\Data\Haixing\Reprocessed\Tumour';
if reProcPath == 0
    return
end
reList = fileFinderAll(reProcPath,'mat');
reList(1,:) = [];

% Now we loop through the files...
for n = 1:numF
    
    % File name of original file...
    orig = list{n,2}

    % Find the best matching file in the other one...
    fcmp = strfind(reList(:,2),orig(1:8));
    fcmp = ~cellfun(@isempty,fcmp);
    fcmp = find(fcmp)
    
    % What if there is no perfect match?
    if numel(fcmp) == 1
        newFile = [reList{fcmp,1} filesep reList{fcmp,2}]
    else
        % Then we pass on this file
        continue;
    end
    
    % Once we have a perfect match, then we can open the files
    tmp = open([list{n,1} filesep list{n,2}]);
    dpn = tmp.dpn;
    
    % Open the other, 'new', file
    tmp = open(newFile);
    new = tmp.dpn;
    
    % Check size
    sz1 = size(dpn.d1.sp)
    sz2 = size(new.d1.sp)
    
    if sz1(1) == sz2(1) && sz1(2) == sz2(2)
        
        dpn.d1.mz = new.d1.mz;
        dpn.d1.sp = new.d1.sp;
        dpn.d1.numPoints = new.d1.numPoints;
        
        % Create a new save name / folder
        savePath = [opPath list{n,2}];
        
        save(savePath,'dpn');
        
    end
        
end

assignin('base','status',status);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dpn,rawFile] = importExistingFile(file)
% Import an H5 or mat file which is to be reprocessed

% Either an h5 or mat file
switch lower(file.ext)

    case 'h5'
        
        % Run the h5 import function
        [mz,x,~,~,~,rawFile] = desiH5LoadFunction(file);
        
        % Place everything into the required dpn structure
        dpn.file    = file;
        dpn.opts    = [];
        dpn.d1.mz   = mz;
        dpn.d1.sp   = x;
        dpn.fig     = [];
        dpn.defP    = dpn.file.dir;
        dpn.mode    = 'single';

    case 'mat'
        
        % Just load the mat file and extract the dpn structure
        tmp = open([file.dir file.nam]);
        dpn = tmp.dpn;
        rawFile.dir = dpn.file.dir;
        rawFile.nam = dpn.file.nam;
        
    otherwise
        error('No otherwise');
end

% Ensure that the rawFile has no extension
dot = strfind(rawFile.nam,'.');
if ~isempty(dot)
    rawFile.nam = rawFile.nam(1:dot(end)-1);
end
rawFile.ext = '';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dpn] = importRawSpectralFile(file,opts)
% Import an imzML/raw file

switch lower(file.ext)
    
    case 'imzml'
        opts = opts.imzml;
        [mz,sp] = desiFunctionIMZML(file,opts);
        
    case 'raw'
        opts = opts.raw;
        [mz,sp,~,~,~,~] = h5waters([file.dir file.nam],opts);
        %mz = [];
        %sp = [];
end

% Now add the parts to the dpn structure
dpn.file    = file;
dpn.opts    = [];
dpn.d1.mz   = mz;
dpn.d1.sp   = sp;
dpn.fig     = [];
dpn.defP    = dpn.file.dir;
dpn.mode    = 'single';
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raw] = determineBestRawFile(file,rawFile,rawType)
% If we are doing reprocessing, then we need to find the best raw/imzML
% file with which to do the reprocessing

% Try to use the one provided
if ~strcmp(rawFile.dir(end),filesep)
    rawFile.dir = [rawFile.dir filesep];
end
p1 = [rawFile.dir rawFile.nam];
p2 = [rawFile.dir rawFile.nam '.' rawType];
if exist(p1) %#ok<EXIST>
    raw = p1;
    return
elseif exist(p2) %#ok<EXIST>
    raw = p2;
    return
end

% Let's try to find the imzML/raw file in the same folder with the
% specified name...
p1 = [file.dir rawFile.nam];
p2 = [file.dir rawFile.nam '.' rawType];
if exist(p1)
    raw = p1;
    return
elseif exist(p2)
    raw = p2;
    return
end

% Now what? Let's try to find any old raw/imzML files in the same folder
list = fileFinderAll(file.dir,rawType);
list(1,:) = [];
if size(list,1) == 1
    raw = [list{1,1} filesep list{1,2}];
    return
elseif size(list,1) > 1
    % We have multiple raw files, so cannot easily decide
    raw = [];
    return
end

% Finally we can look one folder up for a matching file name
[pf] = previousFolder(file.dir);
list = fileFinderAll(pf,rawType);
list(1,:) = [];
if size(list,1) > 0    
    fx = strcmp(list(:,2),file.nam);    
    if sum(fx) == 1
        raw = [pf list{fx,2}];
        return
    end    
end    

% If we get all the way here, then no file matches

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dump
% Use the non-corrected files?
nonCorrect = true;




% Get the options
tmp.ext = lower(rawType);
%[mainOpts,flag] = getOptions(tmp);

if ~flag
    return;
end

% Now we loop through the files...
for n = 1:numF
    
    % Format the file names
    file.dir = [list{n,1} filesep];
    file.nam = list{n,2};
    file.ext = fileExtension(file.nam);
    
    % Add status
    status{n,1} = file.nam;
    
    % Reformat the .dat files to the .raw file name instead
    switch lower(file.ext)
        case 'dat'
            sl = strfind(file.dir,filesep);
            sl = sl(end-1);
            file.nam = file.dir(sl+1:end-1);
            file.dir = file.dir(1:sl);
            file.ext = 'raw';
            
        case 'h5'
            
            % Don't need to do anything!
            
        case 'h5-old'
            sl = strfind(file.dir(1:end-1),filesep);
            sl = sl(end);
            file.nam = file.dir(sl+1:end);
            file.dir = file.dir(1:sl);
            file.ext = 'h5-raw';
            if strcmp(file.nam(end),filesep)
                file.nam = file.nam(1:end-1);
            end
            
            % Here we do something else to ensure that we just get the
            % non-corrected file, from another folder
            if nonCorrect                
                file.dir = [file.dir 'non_corrected' filesep];                
                tmpName = strfind(file.nam,'Analyte 1');
                file.nam = [file.nam(1:tmpName+8) '.raw'];                
            end
    end
        
    % Now decide what to do with it...
    switch lower(file.ext)

        case 'raw'               
            try
                
                % Do the import / processing function here
                [MZ,X,~,~,~,opts] = h5waters([file.dir file.nam],mainOpts.raw);

                % Update the guidata structure and draw a few things here
                dpn.file = file;
                dpn.opts = opts;

                dpn.d1.mz = MZ;
                dpn.d1.sp = X;

                dpn.fig = [];
                dpn.defP = dpn.file.dir;

                dpn.mode = 'single';

                % Now ask the user for a save location...
                saveDir = [dpn.file.dir dpn.file.nam(1:end-6) '-Batch' datestr(now,'yymmdd-HHMMSS') '.mat'];
                save(saveDir,'dpn');
            catch err
                err
                disp(['FAIL WITH FILE: ' file.nam]);
                status{n,2} = 'Import failed';
            end
            
        case 'h5'
            if ~strcmpi(rawType,'imzml')
                error('Do not know what to do');                
            end
            
            % This reads the main parts from the file and thus makes it
            % possible to save everything
            try
                [mz,x,img,opt,allAnno,imzML] = desiH5LoadFunction(file);
                status{n,2} = 'H5 loaded';
            catch
                disp('Skip this file...');
                disp(file.nam);
                status{n,2} = 'H5 loading failed';
                continue;
            end
            
            % We can try to reprocess the data from the imzML file if that
            % is available... Expect to find the imzML folder in the exact
            % location as specified, or somewhere within the 'path' folder
            % as specified by the user'
            poss = false;
            newP = [imzML.dir filesep imzML.nam];
            if exist(newP,'file')
                poss = true;
            else
                newP = [path filesep imzML.nam '.imzML'];
                if exist(newP,'file')
                    imzML.dir = [path filesep];
                    poss = true;                
                end
                
                % Try again
                if ~poss
                    newP = [file.dir file.nam(1:end-3) '.imzML'];
                    if exist(newP,'file')
                        imzML.dir = file.dir;
                        imzML.nam = [file.nam(1:end-3) '.imzML'];
                        poss = true;
                    else
                        newP = [file.dir file.nam(1:end-4) '.imzML'];
                        if exist(newP,'file')
                            imzML.dir = file.dir;
                            imzML.nam = [file.nam(1:end-4) '.imzML'];
                            poss = true;
                        end
                    end
                end
            end
            
            % Potentially reprocess
            if poss
                
                % Expect it to fail
                %try
                    [mz5,x5] = desiFunctionIMZML(imzML,mainOpts.imzml);
                    
                    % Check for sizes...
                    sz = [size(x5,1) size(x5,2)];
                    if size(x,1) == sz(1) && size(x,2) == sz(2)
                        mz = mz5;
                        x  = x5;
                    end
                    status{n,3} = imzML.nam;
                    
                %catch err
                    %err
                    % Shame. It failed
                    %disp('Could not reprocess');
                %end                
            else
                status{n,3} = 'Not reprocessed';
            end
            
            % Now that we have the stuff, let us dump it in the appropriate
            % part and then save the file
            dpn.file = file;
            dpn.opts = mainOpts.imzml;

            dpn.d1.mz = mz;
            dpn.d1.sp = x;
            dpn.d1.img = img;

            dpn.opt.coreg = opt.coreg;

            dpn.fig = [];
            dpn.defP = dpn.file.dir;

            dpn.anno = allAnno;

            dpn.mode = 'single';

            % Temporary re-save location...
            %saveDir = [dpn.file.dir dpn.file.nam filesep f2.nam(1:end-3) 'H2M-NonLockMass.mat'];
            saveDir = [file.dir file.nam(1:end-3) '.mat'];
            save(saveDir,'dpn');
            status{n,4} = [file.nam(1:end-3) '.mat'];
                        
            
        case 'h5-raw'
            
            try
                % Need to get the important bits from the H5 file and put into
                % a format that is suitable for this new toolbox. Though we
                % don't actually want the MZ and X matrix
                f2.dir = [list{n,1} filesep];
                f2.nam = list{n,2};
                [~,~,img,opt,allAnno] = desiH5LoadFunction(f2);

                % Now do the (re)processing
                [MZ,X,~,~,opts] = h5waters([file.dir file.nam],mainOpts.raw);

                % Create all the necessary bits here
                dpn.file = file;
                dpn.opts = opts;

                dpn.d1.mz = MZ;
                dpn.d1.sp = X;
                dpn.d1.img = img;

                dpn.opt.coreg = opt.coreg;

                dpn.fig = [];
                dpn.defP = dpn.file.dir;

                dpn.anno = allAnno;

                dpn.mode = 'single';

                % Temporary re-save location...
                dpn.file.dir = 'H:\Ovarian-NonLockMass\';
                %saveDir = [dpn.file.dir dpn.file.nam filesep f2.nam(1:end-3) 'H2M-NonLockMass.mat'];
                saveDir = [dpn.file.dir f2.nam(1:end-3) '.mat'];
                save(saveDir,'dpn');
            
            catch err
                err
                disp(['FAIL WITH FILE: ' file.nam]);
            end
                
        otherwise
            % Error - there are no other supported filetypes
            disp('File not supported');
    end
    
    clear dpn
end

assignin('base','status',status);
status

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



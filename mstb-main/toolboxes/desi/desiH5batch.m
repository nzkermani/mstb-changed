function desiH5batch
%desiH5batch - convert a folder's worth of h5 files to a mat. This also
%incorporates the metaspace data too, which is saved as part of the main
%file

if ismac
    h5Locn = '/Volumes/Data/Data/Breast/DESI/Sabines Set/';
    %h5Locn = '/Volumes/JSM/DB/LNDatabase/';
    %h5Locn = '/Volumes/JSM/DB/ColorectalFull/';
    %h5Locn= '/Volumes/Data/Data/Ovarian/Luisa/All h5 neg/Batch1_2/';
    %h5Locn= '/Volumes/Data/Data/Ovarian/Luisa/All h5 pos/';
    %saveLocn = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/withHist/';
    saveLocn = '/Users/jmckenzi/Desktop/MTSP/';
    mtsp = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/Engine Dump/MTSP/';
    
else
    h5Locn = 'Z:\Data\Breast\DESI\Sabines Set\';
    saveLocn = 'E:\Dropbox\Imperial\Projects\Metaspace\Engine Dump\withHist\';
    mtsp = 'E:\Dropbox\Imperial\Projects\Metaspace\Engine Dump\MTSP\';
end


mtfile = fileFinderAll(mtsp,'mat',true);

% Find all files in the first folder...
allF = fileFinderAll(h5Locn,'h5',true);
numF = size(allF,1);

% Loop through
for n = 1:numF
        
    % Create as a structure
    file.dir = [allF{n,1} filesep];
    file.nam = allF{n,2};
    file.ext = fileExtension(file.nam);
    
    disp([int2str(n) ' - ' file.nam]);
        
    % New file name?
    nfn = [saveLocn allF{n,2}(1:end-3) '.mat'];
    
    % If it exists, then we skip, as no need to duplicate
    if exist(nfn,'file')
        disp('SKIP');
        continue;
    end
    
    % What are the various potential name systems?
    str2match = {[file.nam(1:end-3) '.mat'];...
        ['ICL--' file.nam(1:end-3) '.mat'];...
        ['ICL--' file.nam(1:end-3) '-centroid.mat']};

    % Match to the various things...
    ffx = false(size(mtfile,1),3);
    ffx(:,1) = strcmp(mtfile(:,2),str2match{1});
    ffx(:,2) = strcmp(mtfile(:,2),str2match{2});
    ffx(:,3) = strcmp(mtfile(:,2),str2match{3});    
    matches = sum(ffx,1) == 1;
    
    % Check that matches works well
    if sum(matches) == 1
        ffx = ffx(:,matches);
        ffx = find(ffx);
    else
        disp('--- no (single) match found');
        continue;
    end        
    
        
    % Import
    try

        % Read from H5
        [mz,x,img,opt,anno,imzML] = desiH5LoadFunction(file);
    
        % Create structure
        dpn.file = file;   
        dpn.opts = 'h5 import';
        dpn.d1.mz = mz;
        dpn.d1.sp = x;
        dpn.d1.img = img;
        dpn.opt.coreg = opt.coreg;
        dpn.defP = saveLocn;
        dpn.anno = anno;
        dpn.mode = 'single';
    
        % Open the correct MTSP file
        mtFile = open([mtfile{ffx,1} filesep mtfile{ffx,2}]);
        
        [dpn2] = mtspMatch(dpn,mtFile);
        dpn.mtsp = dpn2.d1;
        
        % Save
        save(nfn,'dpn');
    
        clear dpn;
        
    catch err
        disp(['FAIL: ' file.nam]);
        err
        disp('---------------------------------------');
    end
        
    
end

end


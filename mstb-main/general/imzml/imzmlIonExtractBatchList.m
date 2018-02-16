function [images] = imzmlIonExtractBatchList(meta,mz,ppm)
%imzmlIonExtractBatch - extract ions from a series of imzML files, but
%using a list for extracting only certain ones and normalising the
%intensities according to global scaling factors determined using the stats
%toolbox...

% Define the folder containing the imzML files
fold    = '/Volumes/JSM/DB/Colorectal DESI/JLA imzML Recalibrated/';
fold2   = '/Volumes/JSM/DB/Colorectal DESI/JLA Reprocess Recalibrated/';

% File names are to be taken from the meta.fileID parameter
numF = size(meta.fileID,1);

% How are we going to store the results?
images = cell(numF,7);

% Sort...
mz = sort(mz);

% Loop through...
for n = 1:numF
    
    disp([int2str(n) '/' int2str(numF) ' - ' meta.fileID{n}]);
    
    % What about doing the same from the mat file of processed data?
    file = [fold2 meta.fileID{n} '.mat'];
    tmp = open(file);
    dpn = tmp.dpn;
    clear tmp;
    
    % Now get the optical image (lo-res) and the ion images
    [dpnImg] = imagesFromDPN(dpn,mz,ppm);
    images{n,6} = dpnImg ./ meta.scaleFac(n);
    
    % H&E image?
    tmpImg = dpn.opt.coreg;
    sz = size(tmpImg);
    if min(sz(1:2)) > 1000
        scFac = 1000 / min(sz(1:2));
        images{n,7} = imresize(tmpImg,scFac);
    else
        images{n,7} = tmpImg;
    end

    % What is the filename?
    file = [fold meta.fileID{n} '-RECAL.imzML'];
    
    % Save the standard information
    images{n,3} = mz;
    images{n,4} = file;    
    
    % Run the extraction function
    [images{n,1},images{n,2},images{n,5}] = extractFromIMZML(file,mz,ppm);    
    images{n,5} = imzmlDate(file);
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp,mz,date] = extractFromIMZML(file,mz,ppm)


[sp,mz] = imzmlIonExtract(file,sort(mz),ppm);    
date = imzmlDate(file);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img] = imagesFromDPN(dpn,mz,ppm)
% Extrac the images as appropriate...

sz = size(dpn.d1.sp);
numI = numel(mz);

img = zeros(sz(1),sz(2),numI);

ppmDev = ppm * mz / 1e6;

% Loop through
for n = 1:numI
    
    % mz values
    mzLo = mz(n) - ppmDev(n);
    mzHi = mz(n) + ppmDev(n);
    
    % Find it...
    mask = dpn.d1.mz >= mzLo & dpn.d1.mz <= mzHi;
    
    % If we find one, then take it...
    if sum(mask) == 1
        % Add to matrix
        img(:,:,n) = dpn.d1.sp(:,:,mask);
    elseif sum(mask) > 1
        % Find the nearest...
        disp('MULTIPLE');
        ff = find(mask);
        gg = mz(ff);
        
        [~,df] = min(abs(gg - mz(n)))
        
        img(:,:,n) = dpn.d1.sp(:,:,ff(df));
    else
        % No image found
    end
    
    
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
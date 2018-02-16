function [images] = desiIonTissueExtract(meta,mz,ppm)
%desiIonTissueExtract - a function that is very complicated. Requires meta
%structure and scaleFac obtained from PQN-Global normalisation option in
%the stats toolbox...

% Define the folder containing the imzML files
fimz = '/Volumes/JSM/DB/Colorectal DESI/JLA imzML Recalibrated/';
fdpn = '/Volumes/JSM/DB/Colorectal DESI/JLA Reprocess Recalibrated/';
probThreshold = 0.95;

% File names are to be taken from the meta.fileID parameter
numF = size(meta.fileID,1);

% How are we going to store the results?
images = cell(numF,10);

% Sort...
mz = sort(mz);

% Loop through...
for n = 1:numF
    
    % Mat file name
    mfn = [fdpn meta.fileID{n} '.mat'];
    
    % Open the mat file
    tmp = open(mfn);
    dpn = tmp.dpn;
    clear tmp;
    
    % Prepare it for segmentation
    [dx,~,anno] = dpnPrepFullMS(dpn,[1 1000],[]);
    [dx] = xxxSegMethodMMC(dx,[],anno);
    
    % Create the temporary image
    sz = size(dpn.d1.sp);
    tmpImg = xxxReturn2Image(real(dx.data),dx.idx,sz(1:2));
    
    % Whcih one of these images should be used as the mask for this tissue
    % type that we are interested in?
    unqC = unique(anno.histID(anno.mask2 > 0));
    cmpH = strcmp(unqC,meta.histID{n});
    maskImg = tmpImg(:,:,cmpH) > probThreshold;
   
    % Now get the optical image (lo-res) and the ion images
    [dpnImg] = imagesFromDPN(dpn,mz,ppm);

    % Also consider extracting from the imzML file as well... Will need to
    % lop the bottom/top row off the imzML images to match the processing
    % performed on them to generate the mat files
    imzFile = [fimz meta.fileID{n} '-RECAL.imzML'];
    [imzImg,imzMZ,imzDate,imzTIC] = extractFromIMZML(imzFile,mz,ppm);
    imzMZ = imzMZ(1:end-1,:,:);
    imzImg = imzImg(1:end-1,:,:);
    imzTIC = imzTIC(1:end-1,:);
    
    % Shall we obtain the H&E image whilst we are here?
    heImg = dpn.opt.coreg;
    sz = size(heImg);
    if min(sz(1:2)) > 1000
        scFac = 1000 / min(sz(1:2));
        heImg = imresize(heImg,scFac);            
    end
    
    % Now let us save the various parts into a cell structure array thing
    % for the visualisation function to use
    images{n,1} = meta.fileID{n};
    images{n,2} = meta.Final_label_for_analysis{n};
    images{n,3} = mz;
    images{n,4} = dpnImg;
    images{n,5} = imzDate;
    images{n,6} = imzImg;
    images{n,7} = imzMZ;
    images{n,8} = imzTIC;
    images{n,9} = maskImg;
    images{n,10}= heImg;
    
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp,mz,date,tot] = extractFromIMZML(file,mz,ppm)

% Extract the images from the imzml file
[sp,mz,tot] = imzmlIonExtract(file,sort(mz),ppm);

% Get the date for the sake of it
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
        [~,df] = min(abs(gg - mz(n)));        
        img(:,:,n) = dpn.d1.sp(:,:,ff(df));
    else
        % No image found
    end
    
    
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
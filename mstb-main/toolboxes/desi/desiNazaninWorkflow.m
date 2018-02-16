function [mz,sp,bin] = desiNazaninWorkflow(path,file,opts)
% desiNazaninWorkflow - function to run Nazanin's version of the code for
% peak alignment and spectral processing. Some of the options are specified
% here as they have not been incorporated into the default options provided
% prior to processing.
%
% INPUTs
% path  - text string of file path
% file  - text string of file name
% opts  - structure of suitable options (to be listed here...)
%       - ppm (number)
%       - binOfLengthOneOK (logical)
%       - splitCorrection (logical)
%
% OUTPUTs
% SPnew     - data matrix
% MZnew     - m/z values
% bin       - details of peak matching (?)
% SizeMZ    - number of mz values (presumably obsolete)
%
% AUTHORs
% Nazanin Z Kermani     | Imperial College London
% James S McKenzie      | Imperial College London
%
% REFs
% a) Kazmi, S. A. et al., Alignment of high resolution mass spectra: 
% Development of a heuristic approach for metabolomics. Metabolomics 2, 75
% 83 (2006).
% b) Kermani et al., (to be published)
%
%
% ISSUEs (various problems and issues with the code)
% - Need to define the options
% - getMSImageProfileMatrix is an incorrect function
% - Need to specify the options for that function, too
% - Log and normalise functions have no effect, as 'data' unused afterwards
% - Solvent peak removal depends on an H5 file, which is stupid as these
% don't exist for unprocessed files. Needs to be entirely re-worked
% - In nzkSolventPeaks, the counter has been replaced, although not tested
% - nzkRandomPeaks requires the background pixels again
% - nzkSpatialRandomness function also full of errors

% Automated definition of the options in here. This needs to be improved
%if nargin == 2
    opts.ppmTol     = 5;
    opts.binLength1 = false;
    opts.splitCorr  = true;
    
    % Spatial randomness options
    opts.spr.minPix = 3;
    opts.spr.numSim = 1000;
    opts.spr.alpha  = 0.05;

%end


% Get the image spectrum - note that this is perhaps the wrong function,
% and could use instead the more traditional one. I think that this
% requires options to be specified as well, which they are not here.
%[mz,sp] = getMSImageProfileMatrix(imzML);
[data] = getimzMLRaw(path,file,opts.mzRange);


return

% Log transformation - what is data. This is badly written code. Data is
% unused after the normalisation stage below...
%data = logData(sp);

% TIC normalisation
%data = ticNorm(data, 'TIC');

% Peak matching
[mz,sp,bin] = nzkPeakMatching(data,opts);

% Exclude solvent peaks
%[mz,sp,bin] = nzkSolventPeaks(mz,sp,bin);

% Spatial randomess
%[mz,sp,bin] = nzkSpatialRandomness(mz,sp,bin,opts)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = getimzMLRaw(path,file,mzRange)
% Simply extract the fully raw profile data from the imzML file and store
% in a cell array

% Add Java package
javaclasspath('packages/imzMLConverter/imzMLConverter.jar')

% Read imzML file
imzML = imzMLConverter.ImzMLHandler.parseimzML([path file]);

% Image dimensions
nCol = imzML.getWidth();
nRow = imzML.getHeight();

data = cell(nRow-1,nCol);
numP = zeros(nRow-1,nCol);

wb = waitbar(0,'imzML');

for y = nRow-1:-1:1
    
    for x = 1:nCol
        
        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end
        
        % Get the data
        mz = imzML.getSpectrum(x,y).getmzArray();
        sp = imzML.getSpectrum(x,y).getIntensityArray();
        
        % Filter the m/z values
        mask = mz >= min(mzRange) & mz <= max(mzRange);
        
        % Save the data
        data{y,x} = [mz(mask) sp(mask)];
        
        numP(y,x) = sum(mask);
    end
    
    waitbar(y/nRow,wb);
end

% Delete the waitbar
delete(wb);

figure; imagesc(numP);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp,bin] = nzkPeakMatching(data,opts)
% Perform peak matching. The output, bin, is a structure containing the
% following elements: mz (m/z ratios) | spectra (the id of the spectrum 
% that mzs belongs to) | centroidspectra (centroid of intensities) |
% centroidmz (centroid of m/zs)

% Run the alignment function
bin = nzkPeakAlignment(data,...
    opts.ppmTol,...
    opts.binLength1,...
    opts.splitCorr);

% Bin to the image, to rebuild the datacube
[mz,sp] = bin2image(mz,sp,bin,'withinImage');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp,bin] = nzkSolventPeaks(mz,sp,bin)
% Exclude solvent peaks. This sub-function is itself unworkable as it
% depends on the H5 file. This is not good data processing. It will have to
% be entirely rewritten, which in turn will produce different results.

% Peaks that are more present in the background
bg_pixels = h5read (strcat([slideFilePath, '.h5']), '/background' );

% number of the tissue object (to) pixels
num_to_pixels = sum(bg_pixels(:)>0);

% number of background pixels
num_bg_pixels = sum(bg_pixels(:)<1);

% Define a counter
%counter = 1;

% Hold the m/z ratios that are proved to be originated by noise
noise = false(size(sp,3),1);
for i = 1:size(sp,3)
    
    a = bg_pixels .* squeeze(sp(:,:,i));
    b = (1-bg_pixels) .* squeeze(sp(:,:,i));
    
    if (sum(sum(length(find(a>0))))/num_to_pixels) <= (sum(sum(length(find(b>0))))/num_bg_pixels)
        noise(i) = true;
        %counter = counter+1;
    end
end

% Exclude noisy ions
sp = sp(:,:,~noise);
mz = mz(~noise);
bin = bin(~noise);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = nzkRandomPeaks(mz,sp,bin,opts)
% Remove random peaks from the spectra

% Exclude ions that occur less than minPixels =3
noise = false(size(sp,3),1);
for i = 1:size(sp,3)
    
    b = (1-bg_pixels) .* squeeze(sp(:,:,i));
    
    if (length(find(b>0)) < opts.minPixels)
        noise(i) = true;        
    end
end

% Exclude the noisy ions
sp = sp(:,:,~noise);
mz = mz(~noise);
bin = bin(~noise);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp,bin] = nzkSpatialRandomness(mz,sp,bin,opts)
% Tests for spatial randomness of ions in the image. First we identify the
% ions that are below the 1% threshold. Then we do the tests themselves

% Define the vector to identify certain ions
spr = false(size(sp,3),1);
% counter = 1;
% noise_ = [];
% numPix = [];
% numTOpixels =  sum(sum( (bg_pixels)));
% precentage = 1;

% Identify the ions - this could be done much more easily by reshaping the
% vector, converting pixels into zero/non-zero etc...
for i = 1:size(sp,3)
    
    b = (bg_pixels) .* squeeze(sp(:,:,i));
    
    if (length(find(b>0)) < (numTOpixels*(precentage/100)))
        noise_(counter) = i;
        numPix(counter) = length(find(b>0));
        counter = counter+1;
    end
end

% Spatial randomness
mz = MZnew;
sp = SPnew;

% spatial randombess test (Clark-Evans test, nearest neighbour criteria)
[index_notNoise, P_value] = clusterCSR2(sp,bg_pixels,...
    noise_,...
    numPix,...
    opts.spr.minPix,...
    opts.spr.numSim,...
    opts.spr.alpha);

% What does this function actually do?
ind_noise = noise_(setdiff(1:length(noise_),index_notNoise));

% Exclude noisy ions
sp = sp(:,:,~noise);
mz = mz(~noise);
bin = bin(~noise);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

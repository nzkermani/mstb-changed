function [SPnew, MZnew, bin, SizeMZ]=workflowDESI(PathName,sample)
%% workflow applys preprocessing on DESI imaging, thermo, centroid mode data 
%%    Input:  PathName: directory of the imzML files
%             sample: name of the imzML file
%            ppm - part per million error of the instrument
%            binOfLengthOneOK - whether bin incluing only one ion is okay
%            (=1) or not (=0)
%    Output:SPnew: ion intencities  
%           MZnew: m/z ratios
%           bin: peak matching out put(more info in the body) 
%           SizeMZ - number of m/z values after preprocessing
%% Author: Nazanin z. Kermani, Imperial College London 2016.
% Reference: 1. Kazmi, S. A.,et al.Alignment of high resolution mass spectra: 
% Development of a heuristic approach for metabolomics. Metabolomics 2, 75–83 (2006).
%            2. to be published 


%% read imzML files
FileName = [sample '.imzML'];
slideFilePath = [PathName, sample];
javaclasspath('packages/imzMLConverter/imzMLConverter.jar')
imzML                = imzMLConverter.ImzMLHandler.parseimzML([PathName, FileName]);
[mzTest,SpTest]      = getMSImageProfileMatrix(imzML);
%% 1. log transformation
% 1.log transform
data =logData(SpTest);
%% 2. Normalization
% 2. TIC normalize
[data] = ticNorm(data, 'TIC');

%% 3. peak matching
% ppm: part per million error
% binOfLengthOneOK : if 1 the bins with only one ion in them stay if 1,excluded 
% splitCorrection : if 1 checks for split bin and correct them, if 0 no split correction is applied;
ppm = 5;
binOfLengthOneOK = 0;
splitCorrection = 1;
% bin is the output of the peak matching algorithm, it is a structure
% containing
%     mz : m/z ratios
%     spectra : the id on the spectrum that m/zs blongs to
%     centroidspectra : centroid of intencities 
%     centroidmz :  centroid of m/zs 
bin = PPMPeakAlignmentv10parB1(mzTest, SpTest, ppm, binOfLengthOneOK, splitCorrection);
% 3.1 bin to image, to rebuild the datacube
[MZnew, SPnew] = bin2image(mzTest, SpTest, bin,  'withinImage');

%% 4. exclude solvent peaks
% Peaks that are more present in the background
bg_pixels = h5read (strcat([slideFilePath, '.h5']), '/background' );
% number of the tissue object (to) pixels
num_to_pixels = sum(bg_pixels(:)>0);
% number of background pixels
num_bg_pixels = sum(bg_pixels(:)<1);

counter = 1;
% holds the m/z ratios that are proved to be originated by noise  
noise_ = [];
for i = 1:size(SPnew,3)
    a = bg_pixels .* squeeze(SPnew(:,:,i));
     b = (1-bg_pixels) .* squeeze(SPnew(:,:,i));
     if (sum(sum(length(find(a>0))))/num_to_pixels) <= (sum(sum(length(find(b>0))))/num_bg_pixels)
        noise_(counter) = i;
        counter = counter+1;
    end
end
% exclude noisy ions
SPnew(:,:,noise_)=[];
MZnew(noise_)=[];
bin(noise_) = [];


% 5. exclude random ions (5.1 happen in few pixels, 5.2 spatialy random)
%% 5.1 exclude ions that occur less than minPixels =3
counter = 1;
noise_ = [];
% minPixels is the minimum number of pixels 
minPixels = 3;
for i = 1:size(SPnew,3)
    b = (1-bg_pixels) .* squeeze(SPnew(:,:,i));
    if (length(find(b>0)) < minPixels)
        noise_(counter) = i;
       counter = counter+1;
    end
end

% exclude noisy ions
SPnew(:,:,noise_)=[];
MZnew(noise_)=[];
bin(noise_) = [];

%% 5.2 find potential noisy spectra (corrected by the test for spatial randomness)
% 5.2.1 Find ions that are present in <1%(precentage) of tissue pixels
counter = 1;
noise_ = [];
numPix = [];
numTOpixels =  sum(sum( (bg_pixels)));
precentage = 1;
for i = 1:size(SPnew,3)
    b = (bg_pixels) .* squeeze(SPnew(:,:,i));
    if (length(find(b>0)) < (numTOpixels*(precentage/100)))
        noise_(counter) = i;
        numPix(counter) = length(find(b>0));
        counter = counter+1;
    end
end

% 5.2.2 Spatial randomness test
mz = MZnew;
sp = SPnew;
% intialize lowestPresence 
lowestPresence = minPixels;
% number of point simulations
numSimulations = 1000;
% level of significane for spatial randomness to retain ions 
levelOfSignificant = 0.05;
% spatial randombess test (Clark-Evans test, nearest neighbour criteria)
[index_notNoise,  P_value] = clusterCSR2(sp,bg_pixels,noise_,numPix,lowestPresence, numSimulations,levelOfSignificant);
ind_noise = noise_(setdiff(1:length(noise_),index_notNoise));
% exclude noisy ions
SPnew(:,:,ind_noise)=[];
MZnew(ind_noise)=[];
bin(ind_noise) = [];

SizeMZ = length(MZnew);

end


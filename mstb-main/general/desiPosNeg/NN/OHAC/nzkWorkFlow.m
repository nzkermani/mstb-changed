function [SPnew, MZnew, bin]=nzkWorkFlow(PathName,sample,ppm,binOfLengthOneOK,...
    splitCorrection,saveSwith,outputPath,method)
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
%%

%% 3. peak matching
% ppm: part per million error
% binOfLengthOneOK : if 1 the bins with only one ion in them stay if 1,excluded
% splitCorrection : if 1 checks for split bin and correct them, if 0 no split correction is applied;
% ppm = 5;
% binOfLengthOneOK = 0;
% splitCorrection = 1;
% bin is the output of the peak matching algorithm, it is a structure
% containing
%     mz : m/z ratios
%     spectra : the id on the spectrum that m/zs blongs to
%     centroidspectra : centroid of intencities
%     centroidmz :  centroid of m/zs
%bin = PPMPeakAlignmentv10parB1(mzTest, SpTest, ppm, binOfLengthOneOK, splitCorrection);
bin = nzkPeakMatching3MethodsV2(path, file, ppm, binOfLengthOneOK,splitCorrection, method);
% 3.1 bin to image, to rebuild the datacube
[MZnew, SPnew] = bin2image(mzTest, SpTest, bin,  'withinImage');

%%
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

if(saveSwith==1)
    save(outputPath,'bin','SPnew','MZnew')
end
end


function [bin,mz,sp] = ohacPeakmatching(path, file, ppm, binOfLengthOneOK,splitCorrection )
%% PPMPeakAlignmentv10par aligns m/z ratios using hierarchical agglomerative clustering
% and ppm based heuristic rules for peak alignement
%%
%    Input:  mz - [rows x columns x mz]
%            sp -  [rows x columns x intensity]
%            ppm - parts per million error of the instrument
%            binOfLengthOneOK - whether bin incluing only one ion is okay
%            (=true) or not (=fals)
%    Output: bin - a structure that for each element of it includes information
%                  of a group ions that represents the same ion species:
%                  mz - mz ratios
%                  Spectra - intensities
%                  centroidmz - centroid mzs
%                  centroidSpectra - centroid intensities
%% Author: Nazanin z. Kermani, Imperial College London 2016.
% Reference: 1. Kazmi, S. A.,et al.Alignment of high resolution mass spectra:
% Development of a heuristic approach for metabolomics. Metabolomics 2, 75–83 (2006).
%            2. Kermani, N. Z.et al, to be published
% Add Java package
javaclasspath('packages/imzMLConverter/imzMLConverter.jar')
fullfile = [path file '.imzML'];
imzML = imzMLConverter.ImzMLHandler.parseimzML(fullfile);

% Get the image spectrum - note that this is perhaps the wrong function,
% and could use instead the more traditional one. I think that this
% requires options to be specified as well, which they are not here.
% Nona What options?
[mz,sp] = getMSImageCentroidMatrix(imzML);
clear imzML
% % Log transformation
sp = sqrt(sp);
% %  normalisation
[sp,~] = jsmNormalise(sp,'tic', 0 , 0);
sclFac = nansum(sp,2);
% replace 0 by 1, nan and inf by 1
sclFac(isinf(sclFac))   = 1;
sclFac(isnan(sclFac))   = 1;
sclFac(sclFac == 0)     = 1;
sp = bsxfun(@rdivide,sp,scaleFac);


%% Build super spectrum (/reference/target)
% Error window is 2 times ppm
Eppm = 2*ppm/1e6;
% pool all the spectra together and sort it out
[sorted_mz_value, sorted_mz_indx] = sort(mz(:));
clear mz
sorted_mz_value = sparse(sorted_mz_value);
indx_zero = min(find(sorted_mz_value>0));
sorted_mz_value(1:indx_zero) = [];
sorted_mz_indx(1:indx_zero) = [];
TempSp = sp(:);
clear sp
TempSp = TempSp(sorted_mz_indx);
% spectra padded with lots of zero's, take the zero's out
% To facilitate paralle alignment
% caculate the distance between consecutive m/z values
diff_mz = diff(sorted_mz_value);

% Scale distances into 2 times ppm error window
eppm = sorted_mz_value*Eppm*2;

% To facilitate paralle alignment
% find gaps in the spectra bigger that 2 times ppm error
% window (to create intervals)
flag_indx = find(diff_mz > eppm(1:end-1));

% creats intervals
clearvars -except TempSp sorted_mz_value  sorted_mz_indx flag_indx ppm binOfLengthOneOK splitCorrection
number_of_intervals = size(flag_indx,1);
flag_indx = [1; flag_indx];
%% Initiate intervals
% Intervals partition the super spectrum and they can be aligned
% independently. Interval is an array of structures that hold the
% information of:
%                mz - m/z rations (ions)
%                I,J,Z - the coordinates on the 3d image cube
%                       (there is scope to delete this (TODO), not important but saves space)
%                spectra - the location on the initial sepectra
i=0;

if(number_of_intervals>0)
    for i=1:(number_of_intervals)
        interval(i).ions = sorted_mz_value((flag_indx(i)+1):(flag_indx(i+1)));
        interval(i).spectra = sorted_mz_indx((flag_indx(i)+1):(flag_indx(i+1)));
        interval(i).intensity = TempSp(((flag_indx(i)+1):(flag_indx(i+1))));
    end
end
interval(i+1).ions = sorted_mz_value((flag_indx(end)+1):end);
interval(i+1).spectra = sorted_mz_indx((flag_indx(end)+1):end);
interval(i+1).intensity = TempSp(((flag_indx(i)+1):end));
clearvars -except number_of_intervals interval ppm binOfLengthOneOK splitCorrection
%% Do the alignment here in parallel for loop

% intiate a parallel pool
%numOfCPUs = getenv('NUMBER_OF_PROCESSORS');
%poolObj = parpool('local', str2double(numOfCPUs));
%poolObj = parpool('local');
tic;
parfor i=1:(number_of_intervals+1)
    % align the intervals individualy
    bin{i} = alignIntervalB1M1(interval(i),ppm, binOfLengthOneOK,splitCorrection) ;
end
toc;
% delete pool object
% delete(poolObj);
%% replace bin structure with bin of type cell to ease the downstream processing
bin = binStruct2bin(bin);
% clear all the variables except bin
clearvars -except bin
end
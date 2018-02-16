function bin = nzkPeakMatching(path, file, ppm, binOfLengthOneOK,splitCorrection)
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
javaclasspath('imzMLConverter/imzMLConverter.jar')
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
sp = bsxfun(@rdivide,sp,sclFac);
dims = size(sp);

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
clearvars -except  dims TempSp sorted_mz_value  sorted_mz_indx flag_indx ppm binOfLengthOneOK splitCorrection
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
clearvars -except  dims number_of_intervals interval ppm binOfLengthOneOK splitCorrection
%% Do the alignment here in parallel for loop

% intiate a parallel pool
%numOfCPUs = getenv('NUMBER_OF_PROCESSORS');
%poolObj = parpool('local', str2double(numOfCPUs));
%poolObj = parpool('local');
tic;
parfor i=1:(number_of_intervals+1)
    % align the intervals individualy
    bin{i} = alignIntervalB1M1(interval(i),ppm, dims, binOfLengthOneOK,splitCorrection) ;
end
toc;
% delete pool object
% delete(poolObj);
%% replace bin structure with bin of type cell to ease the downstream processing
bin = binStruct2bin(bin);
% clear all the variables except bin
clearvars -except bin
end
function bin = nzkPeakAlignment(data,ppm,binOfLengthOneOK,splitCorrection )
% nzkPeakAlignment

% PPMPeakAlignmentv10par aligns m/z ratios using hierarchical agglomerative clustering
% and ppm based heuristic rules for peak alignement
%%
%    Input:  mzTest - [rows x columns x mz]
%            SpTest -  [rows x columns x intensity]
%            ppm - part per million error of the instrument
%            binOfLengthOneOK - whether bin incluing only one ion is okay
%            (=1) or not (=0)
%    Output: bin - a structure that for each element of it includes information
%                  of a group ions that represents the same ion species:
%                  mz - mz ratios
%                  Spectra - intensities
%                  centroidmz - centroid mzs
%                  centroidSpectra - centroid intensities
%% Author: Nazanin z. Kermani, Imperial College London 2016.
% Reference: 1. Kazmi, S. A.,et al.Alignment of high resolution mass spectra:
% Development of a heuristic approach for metabolomics. Metabolomics 2, 75–83 (2006).
%            2. to be published


% Error window is 2 times ppm
Eppm = 2 * ppm / 1e6;

% Need to define all of the mz values into a single vector
sz = size(data);
d2 = reshape(data,[sz(1)*sz(2) 1]);

% Super spectrum of mz values
allSp = vertcat(d2{:});

% Sort all mz values
[sorted_mz_value,sorted_mz_indx] = sort(allSp(:,1));

% Calculate difference between neighbouring m/z values
diff_mz = diff(sorted_mz_value);



% Initiate intervals

% Do the alignment here in parallel for loop

% intiate a parallel pool
%numOfCPUs = getenv('NUMBER_OF_PROCESSORS');
%poolObj = parpool('local', str2double(numOfCPUs));
%         poolObj = parpool('local');

parfor i=1:(number_of_intervals+1)
    %for i=1:(number_of_intervals+1)
    % align the intervals individualy
    bin{i} = alignIntervalB1(interval(i),mzTest,SpTest, ppm, binOfLengthOneOK,splitCorrection) ;
end
% delete pool object
%         delete(poolObj);
%% replace bin structure with bin of type cell to ease the downstream processing
bin = binStruct2bin(bin);

% clear all the variables except bin
%clearvars -except bin

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function determineIntervals(mz)
% THis code is too messy for anyone to make sense of it. Needs to be
% tidied.

% Intervals partition the super spectrum and they can be aligned
% independently. Interval is an array of structures that hold the
% information of:
% mz - m/z rations (ions)
% I,J,Z - the coordinates on the 3d image cube
%   (there is scope to delete this (TODO), not important but saves space)
% spectra - the location on the initial sepectra

% Scale distances into 2 times ppm error window
eppm = mz * Eppm * 2;

% Find gaps in the spectra bigger that 2 times ppm error window, in order
% to create intervals
flagIdx = find(diff_mz > eppm(1:end-1));

% Create the intervals
numInt = size(flagIdx,1);
count = 1;
flagIdx = [0; flagIdx];

i=0;

% This is just a hangup from the old code
indx_zero = 1;

if(numInt>0)
    for i=1:(numInt)
        itv(i).ions = mz((indx_zero+flagIdx(i)):(flagIdx(i+1)+indx_zero-1));
        
        [itv(i).I, itv(i).J,itv(i).Z] = ind2sub(size(mzTest)...
            , sorted_mz_indx((indx_zero+flagIdx(i)):(flagIdx(i+1)+indx_zero-1)));
        itv(i).spectra = sorted_mz_indx((indx_zero+flagIdx(i)):(flagIdx(i+1)+indx_zero-1));
    end
end
itv(i+1).ions = sorted_mz_value((indx_zero+flagIdx(end)):end);

[itv(i+1).I, itv(i+1).J,  itv(i+1).Z] = ind2sub(size(mzTest)...
    , sorted_mz_indx((indx_zero+flagIdx(end)):end));

itv(i+1).spectra = sorted_mz_indx((indx_zero+flagIdx(end)):end);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
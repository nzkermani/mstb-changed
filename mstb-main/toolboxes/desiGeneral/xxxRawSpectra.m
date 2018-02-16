function [ output_args ] = xxxRawSpectra(dpn)
%xxxRawSpectra - function for making images and stuff from raw data
%extracted from Waters raw data

% Determine size of the thing
sz = [size(dpn.d1.sp,1) size(dpn.d1.sp,2)]

% Plot the spectra...
numScan = size(dpn.d1.numPoints);
figure; hold on;
for n = 1:100:numScan
    stem(dpn.d1.numPoints{n}(:,1),dpn.d1.numPoints{n}(:,2));
end


end


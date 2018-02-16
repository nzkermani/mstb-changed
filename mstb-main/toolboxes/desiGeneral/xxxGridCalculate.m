function [ax,ay,bx,by] = xxxGridCalculate(imgA,imgB)
%xxxGridCalculate - do this separately so that we can use it for the H5
%load function

% What is the size of the opt image?
szOp = size(imgA);

% What is the size of the ms image
szMS = size(imgB);

% Determine the ratio between the two
ratio = szOp(1:2) ./ szMS(1:2);

% Linear interpolation across the optical image. Also adjust the gridlines
% to intersect neighbouring pixels, rather than intersecting within a pixel
ax = (1:ratio(2):szOp(2)) - (ratio(2) / 2);
ay = (1:ratio(1):szOp(1)) - (ratio(1) / 2);

% Determine the MS grid
bx = (1:szMS(2))-0.5;
by = (1:szMS(1))-0.5;

end


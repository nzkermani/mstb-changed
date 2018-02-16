function [ peaks ] = rawIMZMLalign(mz,avg)
%rawIMZMLalign - determine which peaks will form part of the cmz vector 
% against which the imzML files will be aligned...

% Normalise
norm = bsxfun(@rdivide,avg,sum(avg,1));

% Average spectrum
med = nanmedian(norm,2);

% Perform peak picking on these
[a,b] = hist(med,1000);
thresh = b(2);
mask = med > thresh;

% Create a structure of picked peaks
peaks.mz = mz;
peaks.sp = med;
peaks.thresh = thresh;
peaks.mask = mask;

% Plot
figure; hold on;
plot(mz,norm);
plot(mz,med,'k');
scatter(mz(mask),med(mask),80,'r','o','filled');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ cmz ] = mzCMZ(data,ppmTol)
%mzCMZ - need a new function to determine the cmz for peak alignment
%purposes. The current version is just a little shit, and the smoothing is
%too variable. Need to bin at a resolution that is suitable and perhaps
%intensity driven.

% High and low mz values
numF = size(data,2);
if numF > 1
    mzv = zeros(numF,2);
    for n = 1:numF
        mzv(n,:) = [min(data(n).mz) max(data(n).mz)];
    end
    flag = false;
    
elseif numF == 1
    numF = size(data,1);
    mzv = zeros(numF,2);
    for n = 1:numF
        mzv(n,:) = [min(data{n}(:,1)) max(data{n}(:,2))];
    end    
    flag = true;
end

% Determine the ppm vector
ppmv = ppmVector(floor(min(mzv(:,1))),ceil(max(mzv(:,2))),ppmTol);

% For each spectrum find the mean intensity
means = cell(1,numF);
for n = 1:numF
    
    means{n} = nanmean(data(n).sp,1);
    
end

% Interpolate...
ip = zeros(numF,numel(ppmv));
for n = 1:numF
    
    mask = data(n).mz > ppmv(1) & data(n).mz < ppmv(end);
    
    [a,~] = comp2vec(data(n).mz(mask),ppmv');
    
    ip(n,a) = means{n}(mask);
    
end

% Plot the vector
%figure; hold on;
%plot(ppmv,ip);

% Let's remove variables that appear in only one of the samples. This is
% potentially dangerous, but worth a try
mask = sum(ip > 0,1) == 1;
ip(:,mask) = 0;



% Mean spectrum
ref = nanmean(ip,1);

%plot(ppmv,-ref,'k');

% Local maxima to determine peak positions...
sl = [ref(2:end) Inf];
sr = [Inf ref(1:end-1)];

% Peak maxima
pp = ref > sl & ref > sr;

% Final cmz vector
cmz = ppmv(pp);


end


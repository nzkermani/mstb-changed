function [ output_args ] = spectrumFWHMdetermine( mz,sp,resol )
%spectrumFWHMdetermine - determine the width of largish peaks over the
%range of a mass spectrum

% Find local maxima over a certain threshold
thr = nanmedian(sp) * 50;

l1 = [sp(2:end) Inf];
l2 = [sp(3:end) Inf Inf];

r1 = [Inf sp(1:end-1)];
r2 = [Inf Inf sp(1:end-2)];


% Local maxima
lm = sp > l1 & ...
    sp > l2 & ...
    sp > r1 & ...
    sp > r2 & ...
    l1 > l2 & ...
    r1 > r2 & ...
    sp > thr;

figure; hold on; 
plot(mz,sp);
stem(mz(lm),sp(lm))

% For each of the peaks with intensity x, determine the span of the peaks
% given a specific resolution...
peaks = find(lm);
numP = numel(peaks);

numR = numel(resol)
fwhms = NaN(numP,numR);

for n = 1:numR
    
    % HM of peaks?
    hm = sp(peaks) / 2;
    
    fwhms(:,n) = hm / resol(n);
    
end

hw = fwhms ./ 2;

l1 = [mz(peaks)'-hw(:,1) mz(peaks)'+hw(:,1)]

l1 = fwhms



stem(mz(peaks)-(fwhms/2),sp(peaks))


end


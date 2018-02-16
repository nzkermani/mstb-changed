function [ val ] = localSNR(mz,sp,locmax)
%localSNR - determine local noise and signal levels...

numV = numel(mz);

% Determine the limits of the regions to take
hw = 50;
st = (1:numV) - hw;
fn = (1:numV) + hw;

% Prevent non-real indices
st(st < 1) = 1;
fn(fn > numV) = numV;

% Calculate the intensity percentiles over the spectrum
val = zeros(1,numV);

% Now calculate...
locmax = find(locmax == 1);
for n = locmax
    %val(1,n) = sp(n) / median(sp(st(n):fn(n)));
    %val(2,n) = sp(n) / mean(sp(st(n):fn(n)));
    %val(3,n) = sp(n) / max(sp(st(n):fn(n)));
    val(1,n) = 1 / std(sp(st(n):fn(n))) / mean(sp(st(n):fn(n)));
end

return

% Use the CV to threshold the spectra
mask = val(1,:) < prctile(val(1,:),5);
%mask = masmooth(mz,double(mask),50,'mean') > 0;

figure; hold on; 
plot(mz,sp,'--k','LineWidth',2); 
%stem(mz,sp);
%plot(mz,val,'LineWidth',2);
scatter(mz(mask),sp(mask));


end


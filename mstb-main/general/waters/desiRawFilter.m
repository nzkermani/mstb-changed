function [ new ] = desiRawFilter( mz,fq,rp )
%desiRawFilter - using a fixed resolving power, determine the frequency of
%specific peaks...

% Determine HWHM for each m/z
hw = (mz / rp) / 2;
ll = mz - hw;
hl = mz + hw;

% Convert to indices (considering that mz is fixed Da width)
df = median(diff(mz));
idx = round((mz-mz(1)) / df) + 1;
ill = ceil((ll-mz(1)) / df) + 1;
ihh = floor((hl-mz(1)) / df) + 1;

ill(ill < 1) = 1;
ihh(ihh > numel(mz)) = numel(mz);

% New spectrum
new = zeros(size(fq));

% Loop through
tic
for n = 1:numel(fq)
    
    % Determine window...
    %fx = mz > (mz(n) - hw(n)) & mz < (mz(n) + hw(n));
    
    new(n) = sum(fq(ill(n):ihh(n)));    
end
toc

% figure; hold on;
% plot(mz,fq,'k','LineWidth',2);
% 
% plot(mz,new,'r','LineWidth',2);

end


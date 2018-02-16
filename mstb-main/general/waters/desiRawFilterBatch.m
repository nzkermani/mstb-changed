function [ flt ] = desiRawFilterBatch(mz,sp)
%desiRawFilterBatch

res = [10000:5000:50000];
numR = numel(res);

flt = zeros(numR,numel(sp));

for n = 1:numR
    
    flt(n,:) = desiRawFilter(mz,sp,res(n));
    
end

flt = bsxfun(@minus,flt,[10000:10000:90000]');

figure;
h = plot(mz,flt);
legend(h,int2str(res'))


end


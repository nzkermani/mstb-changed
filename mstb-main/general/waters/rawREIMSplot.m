function rawREIMSplot(sp)
%rawREIMSplot - plot a file's worth of data

numS = size(sp,1);


totSum = zeros(numS,1);

figure; hold on;

for n = 1:numS
    
    stem(sp{n}(:,1),sp{n}(:,2),...
        'Marker','none');
    
    totSum(n,1) = nansum(sp{n}(:,2));
    
end

figure;
stem(totSum);


end


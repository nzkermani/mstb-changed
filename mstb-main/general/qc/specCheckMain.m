function specCheckMain(mz,sp,grp)
%specCheckMain - this function works on a data matrix with muliple classes,
%by running the function once for each of the classes

[unq,~,ind] = unique(grp);
numG = numel(unq);

mask = mz > 210 & mz < 2000;

for n = 1:numG
    
    % Indices of this class
    fx = ind == n;
    
    x = sp(fx,mask);
    
    % Here perform 2-means clustering
    [idx1,~,~,~] = kmeans(x,2,'replicates',100);
    
    sum(idx1 == 1)
    sum(idx1 == 2)
    
    figure; boxplot(nansum(x,2),idx1);
    
    g1 = nanmean(x(idx1 == 1,:),1);
    g2 = nanmean(x(idx1 == 2,:),1);
    figure; hold on;
    plot(mz(mask),x(idx1 == 1,:),'black');
    plot(mz(mask),x(idx1 == 2,:),'blue');
    plot(mz(mask),-g1,'red');
    plot(mz(mask),-g2,'green');
    
    
    %[inc] = specCheck(sp(fx,mask));
    
    %numel(find(inc == 0))
    
end


end


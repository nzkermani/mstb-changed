function [ output_args ] = incrementalPCA(sp,histID)
%incrementalPCA - add in higher intensity variables, and track cluster
%distance in PCA between two classes

% Let's just normalise the variables at the beginning
sp = bsxfun(@rdivide,sp,nansum(sp,2)) * 1000;
os = nanmedian(sp(sp > 0));

% Determine the mean spectrum
avg = nanmean(sp,1);

% Log transformation
sp = log(sp + os);

% We need to determine how to add in variables
[~,ord] = sort(avg,'ascend');

% Reorder the matrix...
sp2 = sp(:,ord);


% Unique groups...
[unq,~,ind] = unique(histID);
numG = numel(unq);

% Number of variables...
numV = size(sp2,2);
dsts = NaN(numV,numG);
stdv = NaN(numV,numG);

% Loop through...
for n = 1:5:numV
    
    [ll,ss,ee] = pca(sp2(:,1:n),'NumComponents',min([n 2]));
 
    % PC1 centroid...
    for r = 1:numG
        
        dsts(n,r) = mean(ss(ind == r));
        stdv(n,r) = std(ss(ind == r));
        
    end
        
    if dsts(n,1) > dsts(n,2)
        dsts(n,:) = fliplr(dsts(n,:));
        stdv(n,:) = fliplr(dsts(n,:));
    end
        
    %dsts(1:n,:)
    
    
end

xx = 1:numV;
fx = ~isnan(dsts(:,1));

lo = bsxfun(@minus,dsts,stdv);
hi = bsxfun(@plus, dsts,stdv);

figure; hold on;
plot(xx(fx),dsts(fx,:),'LineWidth',10);
plot(xx(fx),lo(fx,:),'b');
plot(xx(fx),hi(fx,:),'r');

end


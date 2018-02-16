function [ include ] = specCheck(x,flag)
%specCheck - determine how similar spectra are. Essentially we are looking
%for bi-modal-like distrubutions. It is important to ensure that pixels
%have the same histological identity!

plotFig = false;

if nargin == 1
    flag = false;
elseif nargin == 2
    flag = true;
end

% Create a histogram of intensities - use this to determine if there is a
% problem with too many low intensities
[hb,ib] = hist(nansum(x,2),5);

if plotFig
    figure; bar(ib,hb);
    
    [a,b,c] = pca(x);
    figure; scatter(b(:,1),b(:,2),80,nansum(x,2),'o','filled');

end

if hb(1) > 3*hb(end)
    % Continue...
else
    % Then not really a distribution issue so use all pixels with a
    % non-zero total sum
    include = sum(x,2) > 0;
    %include = true(size(x,1),1);
    return
end


% 2-means clustering using absolute distances
try
    [idx1,~,~,~] = kmeans(x,2,'replicates',100);
catch
    include = sum(x,2) > 0;
    return
end

% Determine the 'high TIC' group
[grpTic,grpStd] = grpstats(nansum(x,2),idx1,{'mean','std'});

% If we are looking at the background, we want low intensity pixels, not
% the high ones which are more likely to be odd specks of tissue
if flag
    [~,grp] = min(grpTic);
else
    [~,grp] = max(grpTic);
end

% What about mean/prediction intervals at 95%?
[~,cip] = grpstats(nansum(x,2),idx1,{'meanci','predci'});

% Let's use the predicted confidence interval to determine whether new
% observations should be included as it is quite conservative
if flag
    thresh = cip(grp,2);
    
    % Vector of inclusion
    include = nansum(x,2) < thresh;

else
    thresh = max([0 cip(grp,1)]);
    
    % Vector of inclusion
    include = nansum(x,2) > thresh;
end


if plotFig
    figure;
    boxplot(nansum(x,2),include);
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

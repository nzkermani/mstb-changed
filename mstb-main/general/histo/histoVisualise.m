function histoVisualise(img,idx)
%histoVisualise - plot the segmented H&E images
%
% Provide the original full res optical image and a single image of cluster
% membership

% What is the size of the image
sz = size(img);

% How many clusters are there?
unq = unique(idx(:));
unq = unq(unq > 0);
numC = numel(unq);

% How many plots to be made?
numP = numC + 2;
pw = 1/numP;
ph = 0.9;

% All axes
fig = figure('Units','normalized',...
    'Position',[0.03 0.3 0.99 0.4]);
ax = zeros(numP,1);
ax(numC+1,1) = axes('Parent',fig,...
    'Units','normalized',...
    'Position',[pw*(1-1) 0 pw ph]);
imagesc(img);

ax(numC+2,1) = axes('Parent',fig,...
    'Units','normalized',...
    'Position',[pw*(2-1) 0 pw ph]);
imagesc(idx);

% Loop through each cluster
for n = 1:numC
    
    % Indices of pixels of the nth cluster
    fx = idx == n;
    fx = repmat(fx,[1 1 3]);
    
    % Create new image
    img2 = img;
    img2(~fx) = 0;
    
    % Plot it
    ax(n,1) = axes('Parent',fig,...
        'UNits','normalized',...
        'Position',[pw*(2+n-1) 0 pw ph]);
    imagesc(img2);
    
    
end

set(ax,'XTick',[],'YTick',[]);
linkaxes(ax,'xy');

end


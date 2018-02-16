function [ output_args ] = mtspImageObs(data,grp,mz)
%mtspImageOps - make a function that groups spectra in an image etc...

% Best to sort the observations by grp
[~,idx] = sort(grp);
grp = grp(idx);

% If we have specified the m/z value then sort by that too...
if nargin == 3
    [~,varIdx] = sort(mz);
else
    varIdx = 1:size(data,2);
end

% Sort the data matrix accordingly
data = data(idx,varIdx);

% Determine unique
[unq,~,ind] = unique(grp);
unq
fig = figure; 

ax(1) = axes('Parent',fig,...
    'Units','normalized',...
    'Position',[0.05 0.05 0.025 0.9]);
imagesc(ind);
axis off

ax(2) = axes('Parent',fig,...
    'Units','normalized',...
    'Position',[0.1 0.05 0.85 0.9]);
imagesc(data);
axis off

set(ax,'XTick',[],...
    'YTick',[]);

end


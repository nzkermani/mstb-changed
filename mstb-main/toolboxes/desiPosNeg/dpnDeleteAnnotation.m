function dpnDeleteAnnotation(~,~,fig,man)
%dpnDeleteAnnotation - click in the window, determine the nearest
%annotation and then delete it from the window and also the table...

% Get the guidata
dpn = guidata(fig.fig);
if ~isfield(dpn,'anno')
    return
end

% Find all the annotations and then calculate the distance between them all
numP = size(dpn.anno,1);
if numP == 0
    return
end

% Ginput a single click in the axes
axes(dpn.fig.ax.opt(1));
[x,y] = ginput(1);

% Loop through and calculate the distances
allD = zeros(numP,1);
for n = 1:numP
        
    % Determine the centroid of the object
    centx = mean(dpn.anno{n,6});
    centy = mean(dpn.anno{n,7});
    
    % Distance between the click and the object
    allD(n,1) = (centx - x)^2 + (centy - y)^2;
    
end

% Find the nearest...
[~,idx] = min(allD);

% Now delete that patch from the image...
px = dpn.anno{idx,1};
delete(px);

% ...and now from the dpn.anno
tmp = dpn.anno;
tmp(idx,:) = [];
dpn.anno = tmp;

% ... quickly update the guidata
guidata(fig.fig,dpn);

% ...and now from the table
dpnUpdateAnnotationTable(dpn,man);

end


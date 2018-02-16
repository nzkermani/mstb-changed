function [scan2D,resx,resy] = get2Dcoord(x,y)
%

%  round off coordinate values 
x = round(x./0.001); 
y = round(y./0.001);  
% calculate minimum and maximum x & y coordinate values
minx = min(x); maxx = max(x); miny = min(y); maxy = max(y); 
% step for x and y coordinates
stepx  = min(diff(unique(x))); 
ycoord = unique(y); stepy  = min(diff(ycoord));
gridy  = miny:stepy:maxy; gridx = minx:stepx:maxx;
nRows  = length(gridy); nCols = length(gridx);
scan2D = NaN(nRows,nCols);

% loop through all x and y coordinates
iCount = 0;
for iy = gridy
    iCount = iCount + 1;
    if ycoord(iCount) ~= iy;
        continue;
    end
    iScans             = find(y == iy);
    [ix,indcs]         = intersect(gridx,x(iScans));
    scan2D(nRows-iCount+1,indcs) = iScans;   
end
resx = stepx * 0.001;
resy = stepy * 0.001;
return;
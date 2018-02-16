function [x,y] = desiReadRawSize(filename)
% h5waters - convert a Waters RAW file to an H5 file in order to be able to
% batch process a series of the files, e.g. overnight

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

% Default initialisation
[p1,p2] = watersPackages(filename);
numS = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Variable initialisation
xy = zeros(2,numS);
numP = zeros(numS,1);

% Read the scans from the RAW file
for n = 1:numS
    
    % Number of data points in a given spectrum
    numP(n,1) = calllib('MassLynxRaw','getScanSize',p2,1,n);
    
    % Read the xy coordinates for this scan
    xp = libpointer('singlePtr',0);   
    yp = libpointer('singlePtr',0);    
    calllib('MassLynxRaw','getXYCoordinates',p1,1,n,xp,yp);    
    xy(:,n) = [xp.Value yp.Value];
    xy(:,n) = [xp.Value yp.Value];

end

% Determine the xy coordinates for the sample
try
    xy2D = get2Dcoord(xy(1,:),xy(2,:));
catch err
    err
    xy2D = [];
end

% Return the size
[x,y] = size(xy2D);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

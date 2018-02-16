function [scList,numPoints] = watersScans(filename)
% watersScans - determine the number of scans in an image. Hopefully use as
% a quick way to determine the TIC

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

% Draw the waitbar
wb = waitbar(0,'h5waters - initialising');

% Default initialisation
[p1,p2] = watersPackages(filename);

% Number of scans
numS   = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Variable initialisation
xy = single(zeros(2,numS));
totPoints = zeros(numS,2);

% Read the scans from the RAW file
for i = 1:numS
    
    % Number of data points in a given spectrum
    nPoints = calllib('MassLynxRaw','getScanSize',p2,1,i);
    totPoints(i,:) = [nPoints i];
    
    % Read the xy coordinates for this scan
    xp  = libpointer('singlePtr',0);   
    yp  = libpointer('singlePtr',0);    
    calllib('MassLynxRaw','getXYCoordinates',p1,1,i,xp,yp);    
    xy(:,i)  = [xp.Value yp.Value];
    
    % Update the waitbar
    if mod(i,100) == 0
        frac = i/(numS);
        waitbar(frac, wb, ...
            ['RAW > TIC?: ' int2str(i) '/' int2str(numS)],...
            'FontSize',10);
    end
    
end

% Determine the xy coordinates for the sample
[xy2D,~,~] = get2Dcoord(xy(1,:),xy(2,:));

% Size of the image
[nRows,nCols] = size(xy2D);
numPoints = zeros(nRows,nCols);
scanNos = zeros(nRows,nCols);

% Loop through each scan and align the m/z vectors
for x = 1:nRows
    for y = 1:nCols
        if ~isnan(xy2D(x,y))
            
            % Put the totPoints into an image
            numPoints(nRows-x+1,y) = totPoints(xy2D(x,y),1);
            
            % And the scan numbers to make life easier
            scanNos(nRows-x+1,y) = totPoints(xy2D(x,y),2);
            
        end
    end
end

% Delete the waitbar
delete(wb);

% Ask the user to draw a rectangle...
fig = figure; imagesc(numPoints);
rct = getrect;

% Round the rectangle to get integer scan numbers
r1 = [floor(rct(1)) ceil(rct(2))];
r2 = ceil(rct(3:4));
r2 = r1 + r2;

% Check that people haven't selected outside the rectangle
if r1(1) < 1
    r1(1) = 1;
end
if r1(2) < 1
    r1(2) = 1;
end

if r2(1) > nCols
    r2(1) = nCols;
end
if r2(2) > nRows
    r2(2) = nRows;
end

% Define in terms of the mask
mask = zeros(nRows,nCols);
mask(r1(2):r2(2),r1(1):r2(1)) = 1;

% Reshape this
mask = reshape(mask,[nRows*nCols 1]);

% And the scans
scans = reshape(scanNos,[nRows*nCols 1]);

% Extract only the scans that we need
scList = scans(mask == 1);

close(fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

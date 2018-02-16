function [sp,numP,xy,xy2D] = desiReadRaw(filename,recal)
% h5waters - convert a Waters RAW file to an H5 file in order to be able to
% batch process a series of the files, e.g. overnight

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

tic

% Draw the waitbar
wb = waitbar(0,'h5waters - initialising');

% Replace cal function in header?
if recal
    [origRecal] = watersRecal(filename,'');
    if length(origRecal) == 0
        recal = false;
    end
end

% Default initialisation
[p1,p2] = watersPackages(filename);
numS = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Variable initialisation
sp = cell(numS,1);
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

    % Quit if no points
    if numP(n,1) == 0
        disp([int2str(n) ' - no points']);
        continue;
    end
    
    % Preallocate variables
    mz  = single(zeros(1,numP(n,1)));    
    int = single(zeros(1,numP(n,1)));
    
    % Read spectrum
    mzp  = libpointer('singlePtr',mz);    
    intp = libpointer('singlePtr',int);    
    calllib('MassLynxRaw','readSpectrum',p2,1,n,mzp,intp);

    % Save the full spectral data
    sp{n,1} = [mzp.Value; intp.Value];
    
    
    % Update the waitbar
    if mod(n,100) == 0
        frac = n/(numS);
        waitbar(frac, wb, ...
            ['RAW > H5: ' int2str(n) '/' int2str(numS)],...
            'FontSize',10);
    end
    
end

% Restore the original cal function in header (if needed)
if recal
    [~] = watersRecal(filename,origRecal);
end

% Determine the xy coordinates for the sample
try
    xy2D = get2Dcoord(xy(1,:),xy(2,:));
catch err
    err
    xy2D = [];
end

% Delete the waitbar
delete(wb);

toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

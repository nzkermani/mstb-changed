function [cube,spec] = extractWatersRaw(filename,ions,ppmTol)
% extractWatersRaw - extract ion images from a single Waters Raw file.
% Based on Kirill and Keith's original function for reading in a raw file
% for the imaging toolbox, this function uses the basic code in order to
% achieve slightly different aims.

% Kirill Veselkov, Imperial College London, 2015
% Keith Richardson, Waters Corporation, 2015
% James McKenzie, Imperial College London, 2016

% Check that we have the library loaded and ready to go...
[p1,p2] = watersPackages(filename);

% Determine the number of scans
nScans = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Create matrices for the extraction of ion images
numI = numel(ions);
extr = zeros(nScans,numI); % will be reshaped later...
xy = single(zeros(nScans,2));

% Determine the mz tolerances for each of the ions
ppm = ppmTol * ions / 1e6;

% Somewhere to store the average spectrum
mzRes = 0.01;
mzVector = 50:mzRes:1000;
mzSt = mzVector(1);
mzFn = mzVector(end);
spec = zeros(size(mzVector));

% Create the waitbar
wb = waitbar(0,'Extracting');

% Somewhere to draw the spectra during testing
%figure;

tic

% Start to read in the data
iterNum = 0;
for iScan = 1:nScans    
    %tic
    
    % number of data points in a given spectrum
    nPoints = calllib('MassLynxRaw','getScanSize',p2,1,iScan); 
    
    % pre-allocate array of zeros for mz values
    mz = single(zeros(1,nPoints));    
    
    % pre-allocate array of zeros for intensnities
    int = single(zeros(1,nPoints));    
    
    % get a pointer
    mzp = libpointer('singlePtr',mz);  
    
    % get a pointer
    intp = libpointer('singlePtr',int); 
    
    % Read spectrum
    calllib('MassLynxRaw','readSpectrum',p2,1,iScan,mzp,intp);
    try
        xtmp = double(mzp.Value');
        ytmp = double(intp.Value');
    catch
        disp(['Skip scan ' int2str(i)]);
        continue;
    end
    
    %stem(xtmp,ytmp);
    
    % Loop through each of the ions and determine the best peak from that
    % range.  By best we mean most intense, as this methodology should
    % match that enacted for the imzML files.
    for p = 1:numI

        % Find the ions within the appropriate range
        fx = xtmp > ions(p)-ppm(p) & xtmp < ions(p)+ppm(p);
        fy = ytmp .* fx;

        if sum(fy) == 0
            continue;
        end

        % If we find something, then let's add the maximum value to the
        % matrix. This isn't very advanced, as there may be two peaks
        % within this window. Can make it more advanced in the future
        [idx,locmax] = max(fy);

        % It would be useful take an average spectrum
        mask = xtmp >= mzSt & xtmp <= mzFn;
        tmp = round((xtmp(mask)-mzSt) ./ mzRes) + 1;
        spec(tmp) = spec(tmp) + ytmp(mask)';

        % Save to the datacube
        extr(iScan,p) = ytmp(locmax);

    end

    % Now determine the xy location of the data point
    xp = libpointer('singlePtr',0);
    yp = libpointer('singlePtr',0);
    calllib('MassLynxRaw','getXYCoordinates',p1,1,iScan,xp,yp);
    xy(iScan,:)  = [xp.Value yp.Value];
       
    % Change the waitbar every 100
    if mod(iScan,100) == 0
        frac = iScan/nScans;
        waitbar(frac, wb,'Extracting');
        tt = toc;
        %tt/100
    end
    
end

% Close the wait bar
delete(wb);

% Now we need to convert the 2D series of intensities into a meaningful
% image of variables.  Use the xy location of each pixel, but I'm not sure
% that a simple call to the reshape function will work
[crx,cry] = detGrid(xy);
nR = max(crx);
nC = max(cry);
cube = zeros(nR,nC,numI);

% Now loop through and place the things in the right places
for n = 1:nScans    
    cube(crx(n),cry(n),:) = reshape(extr(n,:),[1 1 numI]);    
end

% Trim spec and return the mz vector in it too
mask = spec > 0;
spec = [mzVector(mask)' spec(mask)'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

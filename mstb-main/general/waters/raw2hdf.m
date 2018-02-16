function [h5N,flag] = raw2hdf(rawPath,rawName,h5Path)
% raw2hdf - read in from RAW and export to HDF5

try
    
    % Easily check that this isn't a Mac
    if ismac
        error('DOES NOT WORK ON MAC');
    end
    
    % Create a flag for errors
    flag = [];

    % Check that we have a filesep at the end
    if ~strcmp(rawPath(end),filesep)
        rawPath = [rawPath filesep];
    end
    
    % Provide output path
    if nargin == 2
        h5Path = rawPath;
    end
    
    % Check that we have a filesep at the end of h5Path
    if ~strcmp(h5Path(end),filesep)
        h5Path = [h5Path filesep];
    end
    
    % Generate an H5 file
    dot = strfind(rawName,'.');
    rootN = rawName(1:dot(end)-1);
    h5N = [h5Path rootN '.h5'];
    
    if exist(h5N,'file')
        flag = 'Exists';
        return
    end
    
    % Draw the waitbar
    wb = waitbar(0,'h5waters - initialising');
    
    % Default initialisation
    [p1,p2] = watersPackages([rawPath rawName]);
    numS = calllib('MassLynxRaw','getScansInFunction',p2,1);
    
    % Variable initialisation
    xy = zeros(numS,2);
    totInt = zeros(numS,1);
    totPts = zeros(numS,1);
    
    % Read the scans from the RAW file
    for i = 1:numS
        
        % Number of data points in a given spectrum
        nPoints = calllib('MassLynxRaw','getScanSize',p2,1,i);
        
        % Preallocate variables
        mz  = single(zeros(1,nPoints));
        int = single(zeros(1,nPoints));
        
        % Read spectrum
        mzp  = libpointer('singlePtr',mz);
        intp = libpointer('singlePtr',int);
        calllib('MassLynxRaw','readSpectrum',p2,1,i,mzp,intp);
        
        % Extract...
        mzTmp = mzp.Value;
        spTmp = intp.Value;
        totInt(i,1) = sum(spTmp);
        
        % Remove zero intensities
        fx = spTmp > 0;
        mzTmp = mzTmp(fx);
        spTmp = spTmp(fx);
        
        totPts(i,1) = numel(mzTmp);
        
        % Read the xy coordinates for this scan
        xp  = libpointer('singlePtr',0);
        yp  = libpointer('singlePtr',0);
        calllib('MassLynxRaw','getXYCoordinates',p1,1,i,xp,yp);
        xy(i,:)  = [xp.Value yp.Value];
        
        % Check that there is actually some data to be recorded here and
        % write all of this to an h5 file
        scanName = ['/raw/scan/' int2str(i)];
        
        if totPts(i,1) <= 10
            h5create(h5N,scanName,[max([totPts(i,1) 1]) 2]);
            
        else
            h5create(h5N,scanName,[totPts(i,1) 2],...
                'Deflate',3,...
                'ChunkSize',[10 2]);
        end
        
        % Only write good data, i.e. more than 0 points
        if totPts(i,1) > 0
            h5write(h5N,scanName,[mzTmp' spTmp']);
            h5writeatt(h5N,scanName,'numP',totPts(i,1));
            h5writeatt(h5N,scanName,'xy',xy(i,:));
        else
            h5write(h5N,scanName,[0 0]);
            h5writeatt(h5N,scanName,'numP',0);
            h5writeatt(h5N,scanName,'xy',xy(i,:));
        end
        
        % Update the waitbar
        frac = i/(numS);
        waitbar(frac, wb, ...
            ['RAW > H5: ' int2str(i) '/' int2str(numS)],...
            'FontSize',10);
        
        
        % Clear the memory - though perhaps unnecessary as it gets reused
        clear mzp intp xp yp
    end
    
    % Determine the xy coordinates for the sample
    xy2D = get2Dcoord(xy(:,1)',xy(:,2)');
    [nRows,nCols] = size(xy2D);
    
    h5create(h5N,'/xy2D',[nRows nCols]);
    h5write(h5N,'/xy2D',xy2D);
    
    % Determine totInt / totPts
    totInt2 = zeros(nRows,nCols);
    totPts2 = zeros(nRows,nCols);
    for j = 1:nRows
        for k = 1:nCols
            totInt2(j,k) = totInt(xy2D(j,k));
            totPts2(j,k) = totPts(xy2D(j,k));
        end
    end
    
    h5create(h5N,'/totInt',[nRows nCols]);
    h5write(h5N,'/totInt',totInt2);
    
    h5create(h5N,'/totPts',[nRows nCols]);
    h5write(h5N,'/totPts',totPts2);
    
    % Delete the waitbar
    delete(wb);
    
catch err
    err
    flag = err;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

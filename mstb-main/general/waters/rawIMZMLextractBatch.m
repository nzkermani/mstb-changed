function rawIMZMLextractBatch(files,peaks,ppmTol)
%rawIMZMLextractBatch - use the picked peaks to extract data from the imzML
%files

% Get all files
numF = size(files,1);

% Generate a name for the h5 file
if ismac
    h5Name = ['/Users/jmckenzi/Documents/Box Sync/H5 Proc/Mac-IMZML-Test-' datestr(now,'yymmdd-HHMMSS') '.h5'];
else
    h5Name = ['E:\Box Sync\H5 Proc\Windows-H5-Test-' datestr(now,'yymmdd-HHMMSS') '.h5'];
end

% Define the peaks that were picked using the mask function
peakMZ = [peaks.mz(peaks.mask) peaks.sp(peaks.mask)];
numP = size(peakMZ,1);

% First let's save the mz vector that is common to all of the files
h5create(h5Name,'/mz',size(peakMZ,1));
h5write(h5Name,'/mz',peakMZ(:,1));

% Waitbar
wb = waitbar(0,'File Extraction');

for n = 1:numF
    
    % Create temporary pks structure
    %pks.mz    = list.origMZ(:,n);
    %pks.span  = list.pks.span';
    
    % Name of the file
    fName = [files{n,1} filesep files{n,2}];
    
    % Extract function
    [data,mzvs] = rawIMZMLbatch(fName,peakMZ,ppmTol);
    %data = rand(35,35,numP);
    %mzvs = rand(35,35,numP);
    
    % Need to save the 'data' and other things to an H5 file
    gpname = ['/file' int2str(n)];
    
    % Save the spectral image matrix
    h5create(h5Name,[gpname '/img'],size(data),...
        'Deflate',3,...
        'ChunkSize',[10 10 10]);
    h5write(h5Name,[gpname '/img'],data);
    
    % Save the m/z values for diagnostic purposes
    h5create(h5Name,[gpname '/mzVals'],size(mzvs),...
        'Deflate',3,...
        'ChunkSize',[10 10 10]);
    h5write(h5Name,[gpname '/mzVals'],mzvs);
    
    tic = nansum(data,3);
    h5create(h5Name,[gpname '/tic'],size(tic));
    h5write(h5Name,[gpname '/tic'],tic);
    
    h5writeatt(h5Name,gpname,'filePath',files{n,1});
    h5writeatt(h5Name,gpname,'fileName',files{n,2});
    
    waitbar(n/numF,wb);   

end

delete(wb);


end


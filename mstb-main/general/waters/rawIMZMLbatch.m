function [data,mzvs] = rawIMZMLbatch(file,peaks,ppmTol)
%rawH5mean - determine mean spectrum

wb = waitbar(0,'Hello');

% Gather the options from the structure
[opts] = getOpts(ppmTol);

% Add the package path
%javaclasspath('packages/imzMLConverter/imzMLConverter.jar');

% Get the handle and then size of the imzML file
imzML = imzMLConverter.ImzMLHandler.parseimzML(file);
nCol = imzML.getWidth();
nRow = imzML.getHeight();

% Create an empty image
numP = numel(peaks);
data = zeros(nRow-1,nCol,numP);
mzvs = NaN(size(data));

% Loop through each pixel
for y = 1:nRow-1%nRow-1:-1:1
    for x = 1:nCol
           
        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end
            
        % Get the data
        mz = imzML.getSpectrum(x,y).getmzArray();
        sp = imzML.getSpectrum(x,y).getIntensityArray();
        
        % Here we need to align the spectrum to the peaks/cmz vector
        [j,k] = samplealign2(...
            peaks,...
            [mz sp],...
            'Band',opts.handBand,...
            'Gap',opts.handGap,...
            'Distance',opts.handDist,...
            'Quantile',[],...
            'SHOWCONSTRAINTS',opts.boolSC,...
            'SHOWNETWORK',opts.boolSN,...
            'SHOWALIGNMENT',opts.boolSA);

        % Put into the data matrix
        data(y,x,j) = sp(k)';
        mzvs(y,x,j) = mz(k)';
                
    end
    waitbar(y/nRow,wb);
end

delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getOpts(ppmTol)

opts.ppmTol = ppmTol;
%opts.mzRes = mzRes;
opts.estimationMethod = 'histogram';
opts.display = false;
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
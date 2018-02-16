function [mz,op,opVal] = spCombineJSM(orbData,resolution,minMZ,maxMZ)
% Takes output from imzML reading script 'Data_extract' and combines all
% spectra from one sample. This is a modified version of Emrys' original
% function...
%
% INPUTs
% orbData       - cell array of orbitrap data
% resolution    - resolution of peaks in PPM
% min_mz        - low m/z value, in Da
% max_mz        - high m/z value, in Da
%
% OUTPUTs
% mz            - binned mz vector spaced according to the ppm value
% op            - matrix of intensities for the file

% In this instance we need to compare the traditional binning approach. It
% is all about validation
if numel(resolution) == 2
    ppmRes = max(resolution);
    dalRes = min(resolution);
    validate = true;    
else
    ppmRes = resolution;
    validate = false;
    opVal = [];
end

% The mz range? This is dependent on a ppm value now, rather than a
% Da-derived m/z difference, e.g. 0.001 Da.
mz = ppmVector(minMZ,maxMZ,ppmRes);

% Define spectral parameters etc...
numM = numel(mz);
numS = size(orbData,2);
op = zeros(numS,numM);
scanSums = zeros(numS,3);

% Define the noise regions (perhaps)
%[noise] = defNoise;

% If we are doing the validation stage for the comparison of the two
% binning methods, then this is the bit to get the second set of matrices
if validate
    
    % Vector of uniformly spaced m/z values
    mzDa = minMZ:dalRes:maxMZ;
    
    % Matrix for intensities
    opDa = zeros(numS,numel(mzDa));
end

% Each of the analyses is one sample of x spectra; process each one by one,
% then put them on a combined m/z scale in the final stage...
for n = 1:numS
    
    % Get data for scan n
    sample = orbData{n};
    
    % How many spectra within here?
    numT = size(sample,2);
    
    % Vector for TIC
    tic  = zeros(numT,1);
    
    % Loop through each spectrum to determine the TIC
    for r = 1:numT        
        tic(r,1) = nansum(sample{r}(:,2));        
    end
    
    % Calculate median TIC
    medSp = median(tic);
        
    % Determine 'good' spectra as those with a TIC greater than half the
    % median value
    keep = tic > (0.5 * medSp);
    keep = find(keep == 1);
    numK = numel(keep);
    
    % New matrix to store the 'good' spectra
    tmp = zeros(numK,numM);
    if validate
        tmp2 = zeros(numK,numel(mzDa));
    end
    
    % Create place to store the sums
    sums = zeros(numK,3);
    
    % Now run through these spectra and interpolate to the common mz
    for r = 1:numK
        
        % Get the scan
        x = sample{keep(r)};
        
        % Mask of mz values to reduce the interpolation stage to only what
        % we are interested in
        mask = x(:,1) >= minMZ & x(:,1) <= maxMZ;
        x = x(mask,:);
        
        % Here we perform the interpolation. There are a range of options
        % that could be checked in order to ensure that the thing is
        % optimised
        ii = interp1(x(:,1),x(:,2),mz,'pchip');
        
        % Add into the tmp matrix
        tmp(r,:) = ii';
        
        % Repeat for the validation using the different mz vector
        if validate
            tmp2(r,:) = interp1(x(:,1),x(:,2),mzDa,'pchip');
        end
        
        % Sum intensity ranges
        faMZ = x(:,1) >= 150 & x(:,1) <= 350;
        liMZ = x(:,1) >= 600 & x(:,1) <= 1000;        
        sums(r,:) = [sum(x(faMZ,2)) sum(x(liMZ,2)) sum(x(:,2))];
        
    end
       
    
    % Take the mean of tmp for the actual output. But why? Why not the few
    % most intense spectra rather than the mean which is dragged down by
    % the early / late ones.  Higher intensity is better?
    op(n,:) = nanmedian(tmp,1);
    
    
    % For comparison...
    if validate
        opDa(n,:) = nanmedian(tmp2,1);
    end
    
    % Get the averages of the raw data TICs
    scanSums(n,:) = nanmedian(sums,1);
end

% Remove NaN and Inf values
op(isnan(op)) = 0;
op(isinf(op)) = 0;

if validate
    opDa(isnan(opDa)) = 0;
    opDa(isinf(opDa)) = 0;
    
    opVal.mz = mzDa;
    opVal.op = opDa;
end
    
% Need to compare to the interpolated data...
faMZ = mz >= 150 & mz <= 350;
liMZ = mz >= 600 & mz <= 1000;        

interpSums = [sum(op(:,faMZ),2) sum(op(:,liMZ),2) sum(op,2)];

scanSums

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FTnn] = defNoise

% Define noise regions...
FTnoiseregions=[197.127000000000;350.717000000000;398.801000000000;...
    477.365000000000;533.774000000000;533.802000000000;...
    542.812000000000;687.417000000000;781.647000000000;...
    792.115000000000;1842.69100000000;1909.46700000000;...
    1984.74200000000;1999.56000000000];

% Round the values
FTn=round(FTnoiseregions);

% Increase region by ±5 Da...
FTnn(:,1)=FTn-5;
FTnn(:,2)=FTn+5;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
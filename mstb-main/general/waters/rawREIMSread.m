function [spec,tot] = rawREIMSread(filename)
% rawREIMSread - just import the various data without data processing or
% other things

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

% Draw the waitbar
wb = waitbar(0,'Raw REIMS - GO GO GO');

% Replace cal function in header?
recal = true;
if recal
    [origRecal] = watersRecal(filename,'');
    if length(origRecal) == 0
        recal = false;
    end
end

% Default initialisation
[p1,p2] = watersPackages(filename);

% Number of scans
numS   = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Variable initialisation
%sp = cell(numS,1);
spec = struct('mz',[],'sp',[]);
tot = zeros(numS,1);

% Read the scans from the RAW file
for i = 1:numS
    
    % Number of data points in a given spectrum
    nPoints = calllib('MassLynxRaw','getScanSize',p2,1,i);
    
    % Preallocate variables
    mz      = single(zeros(1,nPoints));    
    int     = single(zeros(1,nPoints));
    
    % Read spectrum
    mzp     = libpointer('singlePtr',mz);    
    intp    = libpointer('singlePtr',int);    
    calllib('MassLynxRaw','readSpectrum',p2,1,i,mzp,intp);
    
    % Save the raw spectra
    %sp{i,1} = [mzp.Value' intp.Value'];
    spec(i).mz = mzp.Value;
    spec(i).sp = intp.Value;
    tot(i,1) = nansum(spec(i).sp);    
    
    % Update waitbar
    frac = i/(numS);
    waitbar(frac, wb, ...
        ['Raw REIMS - ' int2str(i) '/' int2str(numS)],...
        'FontSize',10);
    
end

% Restore the original cal function in header (if needed)
if recal
    [~] = watersRecal(filename,origRecal);
end

% Delete the waitbar
delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

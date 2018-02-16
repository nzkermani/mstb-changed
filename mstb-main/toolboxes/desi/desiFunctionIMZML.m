function [MZ1,X1,timeStamp] = desiFunctionIMZML(file,opts)
%desiIMZML - load an imzML file!

% Now so let's read in the file
disp([file.dir file.nam]);

% Get the acquisition date and time
timeStamp = imzmlDate([file.dir file.nam]);
%timeStamp = now;

if ~strcmpi(file.nam(end-4:end),'imzml')
    file.nam = [file.nam '.imzML'];
end

% Get the handle and then size of the imzML file
imzML = imzMLConverter.ImzMLHandler.parseimzML([file.dir file.nam]);

% If we provide coordinates for the image rather than using the whole one,
% then we need to work on a way to differentiate between the two approaches
if ~isfield(opts,'roi')
    opts.roi = false;
end

if opts.roi
    nCol = opts.roi.xx(2) - opts.roi.xx(1) + 1;
    nRow = opts.roi.yy(2) - opts.roi.yy(1) + 2;
    
    xVals = opts.roi.xx(1):1:opts.roi.xx(2);
    yVals = opts.roi.yy(1):1:opts.roi.yy(2);
    
    
else    
    % These are the full sized coordinates
    nCol = imzML.getWidth();
    nRow = imzML.getHeight();
    
    yVals = 1:nRow;
    xVals = 1:nCol;
    
end

% Create matrices/cells of the suitable size
numPix      = nCol * (nRow-1);
sp.mz       = cell(numPix,1);%nCol,nRow-1);
sp.counts   = cell(numPix,1);%nCol,nRow-1);
sp.mzcounts = zeros(nCol,nRow-1);
sp.X        = zeros(numPix,1);%nCol,nRow-1);
sp.pixel    = zeros(numPix,2);

i = 0;

%return

% Another waitbar
wb = waitbar(0,'Importing imzML');

% Loop through each pixel
for yy = 1:nRow-1%nRow-1:-1:1
    for xx = 1:nCol
        
        % Increase counter
        i = i + 1;
        
        % Determine y and x
        y = yVals(yy);
        x = xVals(xx);
            
        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end 
            
        % Get the data
        xtmp = imzML.getSpectrum(x,y).getmzArray();
        ytmp = imzML.getSpectrum(x,y).getIntensityArray();
                
        % Round off m/z values
        sp.mz{i,1} = round(xtmp ./ opts.mzRes);        

        % Remove duplicated m/z values after rounding off
        [sp.mz{i,1},sp.counts{i,1}] = remDuplmz(sp.mz{i,1},ytmp);
        
        % Save the information...
        sp.mzcounts(xx,yy) = length(sp.counts{i,1});

        % Pixel location
        sp.pixel(i,:) = [xx yy];        
        
        % Count information
        %sp.counts{x,y} = ytmp(:,2);
        sp.X(i,1) = sum(log2(sp.counts{i,1}(sp.mz{i,1} > 500 ) + 20));
         
    end
    waitbar(y/(nRow+5),wb,'Importing imzML');
end

waitbar((y+1)/(nRow+5),wb,'Processing');

% Now we split up the data into two matrices
sp.mzcounts = reshape(sp.mzcounts,[nCol*(nRow-1) 1]);

% Run the m/z alignment function on the sp matrix
[MZ1,X1] = finalProc(sp,nRow,nCol,opts);

delete(wb);

return

% May also wish to combine complementary ion images here as an extra peak
% alignment stage

% Create a guidata structure, and then save the results, and then display
% them in the right boxes...
dpn.file = file;
dpn.opts = opts;

dpn.d1.mz = MZ1;
dpn.d1.sp = X1;

dpn.fig = fig;
dpn.defP = defP;

dpn.mode = 'single';


% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,val] = remDuplmz(mz,val)

diffmz           = diff(mz)';
[~, mzindcs] = find(diffmz>1);
dubmzindcs       = find(diff(mzindcs)>1);

if ~isempty(dubmzindcs)
    counts_temp = val([1 mzindcs+1]);
    for i = dubmzindcs
        counts_temp(i+1) = sum(val(mzindcs(i)+1:mzindcs(i+1)));
    end
    mz  = mz([1 mzindcs+1]);
    val = counts_temp;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MZ,X] = finalProc(sp,nRow,nCol,opts)

% mzcounts is just the number of data points per spectrum
sp.objmzcounts = sp.mzcounts;

% Cumulative sum of vector lengths
mzobjindcs = cumsum([1 sp.objmzcounts']);
 
% Determine the mz values throughout the sample
objectmzs = zeros(1,sum(sp.objmzcounts));
nObjPixels = sum(sp.X > 1);
for i = 1:nObjPixels
    objectmzs(mzobjindcs(i):mzobjindcs(i+1)-1) = sp.mz{i,1}; 
end

% Find all unique values of the mz vector, then the frequency of each
% throughout the sample
MZ = unique(objectmzs); 
n  = hist(objectmzs,MZ); 
clear objectmzs;

% Filter out mz values that don't appear very frequently
MZ = MZ(n > nCol*(nRow-1) * opts.mzFrac); 
X  = zeros(nRow-1,nCol,length(MZ));

%for y = nRow-1:-1:1
%    for x= 1:nCol

for i = 1:size(sp.pixel,1)
        
        % What is the (x,y) location of this pixel?
        x = sp.pixel(i,1);
        y = sp.pixel(i,2);
        
        % This pixel's mz vector
        mz  = sp.mz{i,1};
        if isempty(mz)
            continue;
        end
        
        % Align mz intensities to a common MZ feature vector
        [~,mzindcs,MZindcs] = intersect(mz,MZ); 
        
        % Put into the new matrix
        X(y,x,MZindcs)    = sp.counts{i,1}(mzindcs);

end

% Restore the mz scale to a meaningful set of numbers
MZ = MZ * opts.mzRes;

% Combine peaks that were split
[X,MZ] = combSplitPeaks2(X,MZ,opts.ppmRes);

% Scale median value to 40 for consistency
%medX = median(X(X~=0));
%X = 40 * X ./ medX; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [MZ,X] = dpnImportLegacy(vals,inds,varargin)
% getMSImageDataMatrix exctracts MS image intensity matrix from..
% imzML data object with a given resolution (default 0.001).
%    Input:  imzML  - the MS image data object of a specimen 
%            mztol  - the mass tolerance of measurements, ...
%                      used to merge peaks that are closer than the mass 
%                      tolerance of measurements
%            mzfrac - the fraction of pixels with a given MS value (0.001) 
%    Output: MZ - mz values 
%            X  - output matrix [rows x columns x mz]
% Author: Kirill A. Veselkov, Imperial College London 2012. 

% Get the default parameters for this function, which needs to use TO/BG
% determination for the better determination of the sample's mz vector
[opts] = getVarArgin(varargin);


% Determine the size of the image
nColumns = max(inds(:,2));
nRows    = max(inds(:,3));

% Preallocation of resources
sp.X          = zeros(nColumns,nRows-1);
sp.mzcounts   = zeros(nColumns,nRows-1);
sp.mz         = cell(nColumns,nRows-1);
sp.counts     = cell(nColumns,nRows-1);

% Loop through each pixel
for y = nRows-1:-1:1
    for x = 1:nColumns

        % This is the correct pixel number
        pixelNumber = find(inds(:,2) == x & inds(:,3) == y);

        % Skip empty scans                        
        if isempty(vals{pixelNumber})
            continue; 
        end
        
        % m/z and count information
        sp.mz{x,y} = round(vals{pixelNumber}(:,1) ./ opts.mzRes);
        sp.counts{x,y} = vals{pixelNumber}(:,2);

        % Remove duplicated mz values after rounding off
        [sp.mz{x,y},sp.counts{x,y}] = remDuplmz(sp.mz{x,y}, sp.counts{x,y});

        % Save the information...
        sp.mzcounts(x,y) = length(sp.counts{x,y});
        
        % There are two possible ways to determine the sp.X matrix which we
        % will use for the determination of tissue and bg pixels. THe
        % originally implemented method is based on the TIC, whilst the new
        % method will use morp open/close operations and the determination
        % of the number of objects too...
        switch lower(opts.method)
            case 'tic'
                sp.X(x,y) = sum(log2(sp.counts{x,y}(sp.mz{x,y} > 500 ) + 20));
                
            otherwise % case 'tobg'
                sp.X(x,y) = sum(log2(sp.counts{x,y} + 20));
        end                
    end
end

% Now we may need to determine the to/bg pixel image
switch lower(opts.method)
    case 'tic'
        % Extract object pixels
        sp.X(sp.X > 0) = 1;
        
    otherwise        
        % Determine the TO/BG pixels
        [tobg,~] = getObjectPixels(sp.X,[],[]);
        
        % Apply the morphological operations
        moOpts.imopen  = opts.morpOpen;
        moOpts.imclose = opts.morpClos;
        moOpts.bigoper = [];        
        [tobg2] = doMorpOperators(tobg,moOpts);
        
        % Now use tobg2 to define the TO pixels
        sp.X = tobg2;       
        
        % Draw a graph with the tobg objects
%         f0 = findobj('Tag','tobg');
%         if isempty(f0)
%             f0 = figure('Units','normalized',...
%                 'Position',[0.25 0.25 0.5 0.5],...
%                 'Tag','tobg'); %#ok<NASGU>
%         else
%             figure(f0);
%         end                
%         subplot(1,2,1); imagesc(tobg);
%         subplot(1,2,2); imagesc(tobg2);
end

% Remove non-tissue specific pixels
sp.objectmz                 = reshape(sp.mz,1,(nRows-1)*nColumns);
sp.X                        = reshape(sp.X,1,(nRows-1)*nColumns);
sp.mzcounts                 = reshape(sp.mzcounts,1,(nRows-1)*nColumns);
sp.objectmz(sp.X == 0)      = []; 
sp.objmzcounts              = sp.mzcounts;
sp.objmzcounts(sp.X == 0)   = [];
mzobjindcs                  = cumsum([1 sp.objmzcounts]);
 
objectmzs = zeros(1,sum(sp.objmzcounts));
nObjPixels = sum(sp.X==1);
for i = 1:nObjPixels
    objectmzs(mzobjindcs(i):mzobjindcs(i+1)-1) = sp.objectmz{i}; 
end

% Common mz scale
MZ = unique(objectmzs); 
n  = hist(objectmzs,MZ); 
clear objectmzs;

% Save mz ratios: if the peak is present in at least % of the 
% total sample matrix
MZ = MZ(n > nColumns*(nRows-1) * opts.mzFrac); 
X  = zeros(nRows-1,nColumns,length(MZ));

for y = nRows-1:-1:1
    for x= 1:nColumns
        
        % This pixel's mz vector
        mz  = sp.mz{x,y}; 
        
        % Align mz intensities to a common MZ feature vector
        [~,mzindcs,MZindcs] = intersect(mz,MZ); 
        
        % Put into the new matrix
        X(y,x,MZindcs)    = sp.counts{x,y}(mzindcs);
    end
end

% Restore the mz scale to a meaningful set of numbers
MZ = MZ * opts.mzRes;

% Combine peaks that were split
switch opts.combSplit
    
    case 'fixed'
        [X,MZ] = combSplitPeaks1(X,MZ,opts.mzRes);
    
    case 'ppm'
        [X,MZ] = combSplitPeaks2(X,MZ,opts.ppmRes);
end

% Scale median value to 40 for consistency
medX = median(X(X~=0));
X = 40 * X ./ medX; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getVarArgin(argsin)
% Defaults for the processing parameters...

% List defaults here
opts.mzFrac     = 0.01;
opts.combSplit  = 'fixed';
opts.mzRes      = 0.001;
opts.ppmRes     = 5;

opts.method     = 'tobg';   % or 'tic' for all pixels...
opts.numS       = 4;        % 4 sections on each slide
opts.morpOpen   = 3;        % morphological openings...
opts.morpClos   = 3;        % ...and closings

for i = 1:2:numel(argsin)
    if strcmpi('mzfrac',argsin{i})
        opts.mzFrac = argsin{i+1};
    elseif strcmpi('combsplit',argsin{i})
        opts.combSplit = argsin{i+1};
    elseif strcmpi('mzres',argsin{i})
        opts.mzRes = argsin{i+1};
    elseif strcmpi('ppmres',argsin{i})
        opts.ppmRes = argsin{i+1};
    elseif strcmpi('method',argsin{i})
        opts.method = argsin{i+1};
    elseif strcmpi('nums',argsin{i})
        opts.numS = argsin{i+1};
    elseif strcmp('morpopen',argsin{i})
        opts.morpOpen = argsin{i+1};
    elseif strcmp('morpclos',argsin{i})
        opts.morpClos = argsin{i+1};
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,binInd] = getObjectPixels(X,nbins,binInd)
% getObjectPixels identifies pixels that differentiates tissue object from background 
% Refernce: N. Otsu, "A threshold selection method from gray-level histogram, 
%           IEEE Trans on System Man Cybernetics 9 (1979), no. 1, 62-66.
% X - grayscale image 
% nbins - number of bins for histogram estimation

if isempty(nbins)
    nbins = 20;
end

[h,hvals] = hist(X(:),nbins);

% calculation of the threshold as described by N. Otsu, A threshold 
% selection method from gray-level histogram, 
% IEEE Trans on System Man Cybernetics 9 (1979), no. 1, 62-66.
L       = length(h);
i       = 1:L;
A       = cumsum(h);
B       = cumsum(h.*i);
u       = B ./ A;
tmp     = (A(L) - A);
v       = (B(L) - B) ./ (tmp + (tmp==0));
F       = A .* (A(L) - A) .* (u-v).^2;

% This is the threshold for the determination...
if isempty(binInd)
    [~,binInd] = max(F);
end

X(X <=hvals(binInd) - (hvals(2)-hvals(1))/2) = 0;
X(X > hvals(binInd) - (hvals(2)-hvals(1))/2) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,val] = remDuplmz(mz,val)

diffmz           = diff(mz)';
[mzuniq, mzindcs] = find(diffmz>1);
dubmzindcs       = find(diff(mzindcs)>1);

if ~isempty(dubmzindcs)
    counts_temp = val([1 mzindcs+1]);
    for i = dubmzindcs
        counts_temp(i+1) = sum(val(mzindcs(i)+1:mzindcs(i+1)));
    end
    mz     = mz([1 mzindcs+1]);
    val = counts_temp;
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,MZ] = combSplitPeaks1(X,MZ,msres)
% combine splitted peaks 
diffmz = diff(MZ);
if nargin < 3
    msres = min(diffmz); 
end
msres = msres + msres./10;

[mzuniq,mzindcs]     = find(diffmz>msres);
dubmzindcs           = find(diff(mzindcs)>1);
[nrows,ncols,nvrbls] = size(X); 
X                    = reshape(X,nrows*ncols,nvrbls);
if ~isempty(dubmzindcs)
    X_temp = X(:,mzindcs);
    for i = dubmzindcs
        X_temp(:,i+1) = sum(X(:,mzindcs(i)+1:mzindcs(i+1)),2);
    end
    
end
X  = reshape(X_temp,nrows,ncols,length(mzindcs));
MZ = MZ(mzindcs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,MZ] = combSplitPeaks2(X,MZ,ppmTol)
% Combine splitted peaks - this version should be tailored towards the
% combination of peaks that are within the mz tolerance in PPM values,
% rather than just a fixed value...


diffmz = diff(MZ);
% if nargin < 3
%     msres = min(diffmz); 
% end
% msres = msres + msres./10;

% Calculate the 5ppm difference
ppm = ppmTol * MZ / 1e6;

% Convert the msres values into a ppm value, which will necessarily be a
% vector for each individual mz value

% Find peaks that are separated from their neighbour by (at least) the 
% ppm difference
[~,mzindcs] = find(diffmz > ppm(1:end-1));

% Determine the difference between these peaks, i.e. those that are
% different than more than the ppm tolerance
diffIND = diff(mzindcs);

% Now find differences between this vector, i.e. references to mzindcs
% which are potentially split
dubmzindcs = find(diffIND > 1);

% If there are no potentially split peaks, then we return having made no
% changes to X and MZ
if isempty(dubmzindcs)
    return;
end

% Size of X
[nrows,ncols,nvrbls] = size(X); 

% Reshape X
X = reshape(X,[nrows*ncols nvrbls]);

% Gather the variables that were determined to be unique and unsplit into
% a new variable.  Do the same for the MZ vector, which can also be trimmed
X_temp = X(:,mzindcs);

% What about taking an alternative MZ vector, that uses the mean of peaks
% that have been split and then re-joined? There is scope for this...
MZ_alt = MZ(mzindcs);
MZ_wgt = MZ(mzindcs);

% Non-complementary list of variables
noncomp = zeros(numel(dubmzindcs),1);

% Loop through each split peak and assign to its neighbouring friend...
% Does this always assign the peak to the best neighbour? It seems to
% always add to the neighbour of higher m/z value?
for i = dubmzindcs
        
    % Check that the intensity patterns are complementary, i.e. they are
    % like merged barcodes with no
    % barc = X(:,mzindcs(i)+1:mzindcs(i+1)) > 0;
    % barc2 = sum(barc,2);    
    % if any(barc2 > 1)
    %     figure; imagesc(barc);
    %     figure; stem(sum(barc,2))
    %     noncomp(i,1) = 1;
    % end

    % This combines peaks that have been split into the new X intensity
    % matrix
    X_temp(:,i+1) = sum(X(:,mzindcs(i)+1:mzindcs(i+1)),2);
    
    % This is a way to calculate a more representative MZ vector.
    % Originally the first variable was used to define the MZ value;
    % however, the combination of multiple later variables is unaccounted
    % for by this method.  By taking the average of the MZ values of the
    % split peaks, a better approximation of the MZ value can be
    % determined. This can be accomplished in two ways:
    %
    % 1) simple mean of MZ values
    MZ_alt(i+1) = mean(MZ(:,mzindcs(i)+1:mzindcs(i+1)));
    %
    % 2) weighted mean according to sparsity of mz values
    pop = X(:,mzindcs(i)+1:mzindcs(i+1)) > 0;
    pop = sum(pop,1);    
    MZ_wgt(i+1) = sum(MZ(:,mzindcs(i)+1:mzindcs(i+1)) .* pop) / sum(pop);
    
    
end

% Finally trim the mz vector
MZ = MZ(mzindcs);

% Which MZ vector to use?
MZ = MZ_wgt;

% Final return X to image size
X  = reshape(X_temp,nrows,ncols,length(mzindcs));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
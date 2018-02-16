function [MZ,X,xy,xy2D,opts] = h5watersV2(filename,opts)
% h5waters - convert a Waters RAW file to an H5 file in order to be able to
% batch process a series of the files, e.g. overnight

% Easily check that this isn't a Mac
if ismac
    error('DOES NOT WORK ON MAC');
end

tic

% Default options, or user-definable options to be included later
if nargin == 1
    [opts] = getDefaults;
end

% Draw the waitbar
wb = waitbar(0,'h5waters - initialising');

% Default initialisation
[p1,p2] = watersPackages(filename);

% Number of scans
numS   = calllib('MassLynxRaw','getScansInFunction',p2,1);

% Variable initialisation
sp = cell(1,numS); 
xy = single(zeros(2,numS));
spcounts = zeros(1,numS);
irobstd  = zeros(1,numS);  
imedian = zeros(1,numS);

tic

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
    
    % Trim out the variables which lie outside the specified m/z range
    mask = mzp.Value' > min(opts.mzRange) & mzp.Value' < max(opts.mzRange);
    mzp.Value = mzp.Value(mask);
    intp.Value = intp.Value(mask);
    
    
    % Here we do peak detection in instances where there has been no lock
    % mass correction performed.
    if opts.peakDetect == 1
        
        % Extract the vectors
        mz1 = double(mzp.Value');
        sp1 = double(intp.Value');
        
        % Resample the spectrum at a determined resolution
        [a,b] = msresample(mz1,sp1,20000,'SHOWPLOT',false);
        b(b < 0) = 0;
        
        % Smooth it
        [bs] = movingWindow(b,5,'mean');        
        bs(bs < 0) = 0;
        
        % Find peaks
        sl = [bs(2:end); Inf];
        sr = [Inf; bs(1:end-1)];
        fx = (b > bs*3) & (b >= sl*3) & (b >= sr*3);
        
        % Draw a plot to show us what we have got to so far
%         figure; hold on;
%         plot(mz1,sp1,'k');
%         plot(a,b,'r');
%         plot(a,bs*3,'b');
%         stem(a(fx),b(fx),'m');
        
        % Now with the rough peak intensities, we need to find the absolute
        % maximum intensity and associated m/z value. This way we preserve
        % the original resolution of the data, even if it is low. Also,
        % this means that as Waters m/z values are constant across scans,
        % we shouldn't need any alignment of m/z values.
        
        allP = find(fx);
        c = zeros(numel(allP),2);
        for r = 1:numel(allP)
            s = allP(r);
            rn = [a(s-1) a(s+1)];
            mask = mz1 >= rn(1) & mz1 <= rn(2);
            [f,g] = max(sp1(mask));
            if ~isempty(f)
                g = g + find(mask,1,'first') - 1;            
                c(r,:) = [mz1(g) f];
            else
                % Do nothing
            end
        end   
        
        % Ditch empty peaks
        cx = c(:,1) ~= 0;
        c = c(cx,:);
           
        %stem(c(:,1),c(:,2),'c');
        
        % Save in the correct format
        sp{i} = c';
    
    else
        % Save the full spectral data
        sp{i} = [mzp.Value; intp.Value];
    end
    
    % Determine mean and std of intensity values
    %irobstd(i)  = 1.4826*median(abs(intp.Value - median(intp.Value)));    
    %imedian(i)  = median(intp.Value);    
    
    % Select only non-zero mz and intensity values
    %nonzeroIndcs = mzp.Value > 0 & intp.Value > 0; %(irobstd(i)/2);
    %sp{i} = sp{i}(:,nonzeroIndcs);
    
    % Round mz values according to the instrument tolerance
    %sp{i}(1,:) = round(sp{i}(1,:) ./ opts.mztol); 
    
    % Remove duplicate mz values
    %[sp{i}(1,:),sp{i}(2,:)] = remDuplmz( sp{i}(1,:), sp{i}(2,:));
    
    % Count the length, i.e. number of variables in this scan
    spcounts(i) = length(sp{i}(1,:));
    
    % Read the xy coordinates for this scan
    xp  = libpointer('singlePtr',0);   
    yp  = libpointer('singlePtr',0);    
    calllib('MassLynxRaw','getXYCoordinates',p1,1,i,xp,yp);    
    xy(:,i)  = [xp.Value yp.Value];
    
    % Update the waitbar
    if mod(i,100) == 0
        frac = i/(numS);
        waitbar(frac, wb, ...
            ['RAW > H5: ' int2str(i) '/' int2str(numS)],...
            'FontSize',10);
    end
    
    % Clear the memory - though perhaps unnecessary as it gets reused
    clear mzp intp xp yp 
end

toc

% Now we need to place into a matrix, which is formed from all values...
all = horzcat(sp{:});
allMZ = unique(all(1,:));

% New matrix
spmat = sparse(numS,numel(allMZ));

% For each scan whack it in the sparse matrix
for n = 1:numS
    
    tmp = unique(sp{n}','rows');
    
    [~,ia,ib] = intersect(allMZ',tmp(:,1));
    spmat(n,ia) = tmp(:,2);
end
    

% Threshold to filter out noisy peaks
peakThr = mean(imedian) + mean(irobstd)*5;

% Start to collate all the m/z values throughout the sample
objectmzs  = zeros(1,sum(spcounts));
mzobjindcs = cumsum([1 spcounts]);

for i = 1:numS
    
    % This is the mz vector for this pixel
    mz = sp{i}(1,:);
    
    % Set intensities less than the threshold to 0
    mz(sp{i}(2,:) <= peakThr) = 0;
    
    % Put the m/z values into a single long vector
    objectmzs(mzobjindcs(i):mzobjindcs(i+1)-1) = mz;
    
    % Save only those non-zero values
    sp{i} = sp{i}(1:2,mz~=0);
    
    % Normalize profiles - why do we have to do this?
    sp{i}(2,:) = sp{i}(2,:)./imedian(i);
end

% Filter out the values with a zero intensity
objectmzs = objectmzs(objectmzs~=0);

% Determine unique m/z values, against which all pixels can be aligned
MZ = unique(objectmzs);

% See how frequently each m/z value occurs
n  = hist(objectmzs,MZ);
clear objectmzs;

% Remove m/z values that appear in too few of the pixels
% Is this done too early?
MZ = MZ(n > numS*opts.mzfrac); 

% Determine the xy coordinates for the sample
xy2D          = get2Dcoord(xy(1,:),xy(2,:));
[nRows,nCols] = size(xy2D);

% Preallocate the spectral datacube. Can we make this sparse without % affecting the rest of the code?
X = zeros(nRows,nCols,length(MZ));

% Loop through each scan and align the m/z vectors
iterNum = numS;
for x = 1:nRows
    for y = 1:nCols
        if ~isnan(xy2D(x,y))
            
            % Extract this scan's m/z vector
            mz  = sp{xy2D(x,y)}(1,:);
            
            % Align mz values to the master MZ list
            [mzjoint,mzindcs] = intersect(mz,MZ);
            [~,MZindcs] = intersect(MZ,mzjoint);
            
            % Place the data in the cube
            X(x,y,MZindcs) = sp{xy2D(x,y)}(2,mzindcs);
            
            % Update the waitbar again...
            if nargin > 4 && mod(iterNum,100) == 0
                frac = iterNum/(2*numS);
                waitbar(frac, wb, ...
                    ['Aligning spectra: ' int2str(iterNum) '/' int2str(2*numS)]);
            end
                    iterNum = iterNum + 1;
        end
    end
end

% Convert the MZ vector to the proper scale
MZ = MZ * opts.mztol;

%assignin('base','Xlow',X);
%assignin('base','MZlow',MZ);

% Combine peaks that are split, i.e. slight mis-alignment problems. This
% could perhaps benefit from a ppm-based approach
[X,MZ] = watersCombSplitPeaks(X,MZ,opts.mztol);

% Perhaps here we run the analysis again but using the image
% complementarity function
[MZ,X] = combCompNeighbours(MZ,X,opts.ppmCombine);

% Determine the overall median intensity value
%medX = median(X(X~=0));

% Scale median value to 40 for consistency - this is not particularly
% necessary I think
%X = 40 * X ./ medX;

% Now that we have the processed data, we should save it to a file...
%newName = [filename(1:end-3) 'mat'];
%save(newName,'MZ','X','xy','xy2D','filename','opts');
%disp(newName);

% Delete the waitbar
delete(wb);

toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getDefaults

opts.mztol      = 0.01;
opts.mzfrac     = 0.01;
opts.peakDetect = 1;
opts.ppmCombine = 10;
opts.mzRange    = [100 1000];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
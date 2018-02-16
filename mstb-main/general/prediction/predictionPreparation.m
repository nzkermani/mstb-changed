function [data] = predictionPreparation(data,opts)
%predictionPreparation - prepare the 'rat' structure to get a set of
%aligned dataset, and chck that the data does not need to be normalised
%
% James McKenzie, 2016

% Let's set the background pixels to NaN values and totally ignore them
[data] = determineTOBG(data,opts);

% Reshape the matrices
[data] = shIm2Mat(data);

% An additional stage in here to determine if the peaks in the files are
% all worth keeping, or if we should be doing some kind of filtration in
% order to remove the peaks which have split / decayed. These will have
% some kind of an impact, but don't know how much just yet...
mzFiltration(data,opts);

%return

% First we need to do the peak matching between the m/z vectors
[data] = mzAlignment(data,opts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = determineTOBG(data,opts)
% Decide if background is to be analysed / predicted

for n = 1:size(data,2)
    
    if opts.remBG
        data(n).tobg = pixTOBG(nansum(data(n).sp,3),[],true);
    else
        data(n).tobg = true(size(data(n).sp,1),size(data(n).sp,2));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = shIm2Mat(x)
% Shape from images to matrices

% The data may already have been reshapen if previously performed
sz1 = size(x(1).sp);
if numel(sz1) == 2
    return;
end

for n = 1:size(x,2)
    
    % Reshape data
    x(n).sz = size(x(n).sp);
    x(n).sp = reshape(x(n).sp,[x(n).sz(1)*x(n).sz(2) x(n).sz(3)]);

    if ~isempty(x(n).anno)
        x(n).anno = reshape(x(n).anno,[x(n).sz(1)*x(n).sz(2) size(x(n).anno,3)]);
    end
    
    if isfield(x(n),'tobg')
        x(n).tobg = reshape(x(n).tobg,[x(n).sz(1)*x(n).sz(2) 1]);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mzFiltration(data,opts)
% Decide which peaks should be kept and which should not.  This is
% necessary perhaps when there is a large discrepancy 

figure; hold on;

% Plot the average specturum for each of the files
numF = size(data,2)
h = zeros(numF,1);
cols = jet(numF);
for n = 1:numF
    
    mn = nanmean(data(n).sp,1);
    mean(mn(:))
    
    [t1,t2] = insertZeros(data(n).mz,mn,0.01,false);
    
    h(n,1) = plot(t1,t2,'-o','Color',cols(n,:));
    
end
    
legend(h);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = mzAlignment(data,opts)
% First define the cmz vector, and then align everything to it

% Define the peak matching options
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

% Prep the mz vectors
numF = size(data,2);
mzPeaks = cell(numF,1);
for n = 1:numF    
    mzPeaks{n,1} = [data(n).mz' mean(data(n).sp,1)'];
end

% Determine the cmz vector
% cmz = mspmatch(mzPeaks,...
%     'estimationMethod','histogram',...
%     'mzRes',0.01,...
%     'display',false);
cmz = mzPeaks{1};

% Now run the samplealign2 function to match A to B...
for n = 1:numF
    
    % Align file n to CMZ
    [j,k] = samplealign2(...
        mzPeaks{1},...
        mzPeaks{n,1},...
        'Band',opts.handBand,...
        'Gap',opts.handGap,...
        'Distance',opts.handDist,...
        'Quantile',[],...
        'SHOWCONSTRAINTS',opts.boolSC,...
        'SHOWNETWORK',opts.boolSN,...
        'SHOWALIGNMENT',opts.boolSA);

    % Make the matrices to finish...
    newSP = zeros(size(data(n).sp,1),size(cmz,1));
    newSP(:,j) = data(n).sp(:,k);
    
    % Save the new matrix and new mz vector
    data(n).mzAl = cmz(:,1);
    data(n).spAl = newSP;
    data(n).szAl = data(n).sz;
    data(n).szAl(3) = size(cmz,1);
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

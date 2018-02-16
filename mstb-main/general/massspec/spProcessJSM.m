function [mz2,var2,opts] = spProcessJSM(mz,op,varargin)
%spProcessJSM - process the variables produced...
%
% Performs a series of operations, such as peak picking to get the data in
% a manageable output format, all hopefully aligned...

% Get the options... Worth reading...
[opts] = getArgsIn(varargin);

% Do RSPA (or not?)
if opts.rspa.rspa
    op = rspa(op,mz,nanmean(op,1),...
        opts.rspa.maxPeakShift,...
        opts.rspa.recursion,...
        opts.rspa.minSegWidth);
end

% Take reference spectrum and smooth it
refSp = nanmedian(op,1);
refSm = masmooth(mz,refSp,opts.smSpan,opts.smType);

% Perform peak picking
defThresh = median(refSm(refSm > 0)) * opts.peakThresh;
pks = msPeakFinder(mz,refSm,'threshold',defThresh);

% Extract new set of variables
[mz2,var2] = varExtract(pks,mz,op,opts.varExtract);

% Normalise
[var2] = jsmNormalise(var2,opts.normMethod,opts.normThresh,0);

% Perform MVA
os = nanmedian(var2(var2 > 0));
[pca.l,pca.s,pca.e] = princomp(log2(var2+os),'econ');

% PCA scatter plots
figure;

subplot(1,2,1); 
scatter(pca.s(:,1),pca.s(:,2),60,1:size(var2,1),'o','filled');
title(['PCA: ' opts.normMethod]);

subplot(1,2,2);
stem(mz2,pca.l(:,1:2));
title('PCA: Loadings');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opts = getArgsIn(argsin)
% Get the options

% Specify defaults here
opts.rspa.rspa = true;         % Decides if you do it or not...
opts.rspa.maxPeakShift = 0.05;  % Da
opts.rspa.recursion = 1;        % True
opts.rspa.minSegWidth = 1;      % (Default)

opts.smType = 'mean';           % Mean, median, etc...
opts.smSpan = 20;               % Smoothing span

opts.varExtract = 'std';        % Use 'sum'/'max', or 'std'

opts.peakThresh = 0.1;       	% Multiplier of the mean, i.e. mean/10

opts.normMethod = 'pqn-mean';   % 'pqn-mean','pqn-median','tic','none'
opts.normThresh = 0;            % Normalisation threshold

nArgs = length(argsin);
for i=1:2:nArgs    
    if strcmpi('rspa',argsin{i})        
        opts.rspa.rspa = argsin{i+1};            
    elseif strcmpi('maxpeakshift',argsin{i})        
        opts.rspa.maxPeakShift = argsin{i+1};            
    elseif strcmpi('recursion',argsin{i})
        opts.rspa.recursion = argsin{i+1};
    elseif strcmpi('minsegwidth',argsin{i})
        opts.rspa.minSegWidth = argsin{i+1};
    elseif strcmpi('smtype',argsin{i})
        opts.smType = argsin{i+1};
    elseif strcmpi('smspan',argsin{i})
        opts.smSpan = argsin{i+1};
    elseif strcmpi('peakthresh',argsin{i})
        opts.peakThresh = argsin{i+1};
    elseif strcmpi('normmethod',argsin{i})
        opts.normMethod = argsin{i+1};
    elseif strcmpi('normthresh',argsin{i})
        opts.normThresh = argsin{i+1};    
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,var] = varExtract(pks,mz,op,method)
% Extract intensities from locations found in pks...

% For each of the mz values specified in pks, we need to find the
% best-matching index in mz
numP = size(pks{1},1);
allP = zeros(numP,1);

% Loop through
for n = 1:numP    
    [~,allP(n,1)] = min(abs(mz - pks{1}(n,1)));    
end

% Allocate space for the new variables
var = zeros(size(op,1),numP);

% Two possible ways to do this:
switch method
    case 'max'
        % This takes the largest intensity one datapoint to either side 
        % of the peak centre. This accounts for flat-tops and peaks that 
        % are a little bit shifted
        for n = 1:numP    
            var(:,n) = nanmax(op(:,allP(n,1)-1:allP(n,1)+1),[],2);
        end
        
    case 'sum'
        % This sums the datapoint and its two nearest neighbours
        for n = 1:numP    
            var(:,n) = nansum(op(:,allP(n,1)-1:allP(n,1)+1),2);
        end
        
    otherwise
        % This just takes the intensity at the exact location
        var = op(:,allP);
end
    
% Extract the mz vector
mz = pks{1}(:,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
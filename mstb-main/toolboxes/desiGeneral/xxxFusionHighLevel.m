function [mzF,spF,hl] = xxxFusionHighLevel(dx,method,compValue)
%xxxFusionHighLevel - give it two data sets and it will do PCA as a
%dimension reduction method and then concatenate the appropriate number of
%PCs and output this new data reduced data set.
%
% James McKenzie, 2016.
%
% Needs to be compatible with QDS and Seg
%
% Note that the data need to have been norm/etc before
%
% dx        - structure containing mz,sp,(idx) for each dataset
% method    - which of the fusion methods [HL (PCA #) | HL (PCA %)]
% compValue - either # or % value

% How many datasets are there? Typically this will be two, but we may wish
% to extend this function to other applications
numD = max(size(dx));

% Determine the maximum number of components to be calculated
switch lower(method(4:end))    
    case '(pca #)'
        maxComp = compValue;
        
    case '(pca %)'
        allSz = zeros(numD,2);
        for n = 1:numD
            allSz(n,:) = size(dx(n).sp);
        end
        maxComp = min(allSz(:));
        
    otherwise
        error('There is no otherwise');
end

% As the data have been norm/etc, we run PCA calculating a limited number
% of components...
hl(numD) = struct('ll',[],'ss',[],'ee',[]);
for n = 1:numD
    
    % This is the raw PCA calculation
    [hl(n).ll,hl(n).ss,hl(n).ee] = pca(dx(n).sp,...
        'NumComponents',maxComp);
    
    % Convert the eigenvalues
    hl(n).ee = 100 * hl(n).ee / sum(hl(n).ee);

    % Scale the variance of the components... This prevents one
    % mode/dataset from dominating the other one
    hl(n).vr = std(hl(n).ss(:,1),[],1);
    hl(n).ss = hl(n).ss / hl(n).vr;
    hl(n).ll = hl(n).ll / hl(n).vr;

end

% Here is where we need to decide on how many components to include
% if doing it based on the % of variance
if strcmpi(method(4:end),'(pca %)')
    
    % What is the maximum value of variance that we want to include?
    maxPerc = compValue;
    numComp = zeros(numD,1);
    for n = 1:numD          
        
        % Cumsum the % variance, find the first over the limit
        numComp(n,1) = find(cumsum(hl(n).ee) > maxPerc,1,'first');
        
        % Check that we have at least two
        if numComp(n,1) == 1
            numComp(n,1) = 2;
        end
    end
    
    % Trim out the crap
    for n = 1:numD        
        hl(n).ll = hl(n).ll(:,1:numComp(n,1));
        hl(n).ss = hl(n).ss(:,1:numComp(n,1));
        hl(n).ee = hl(n).ee(1:numComp(n,1));
    end
    
else
    numComp = repmat(maxComp,[numD 1]);
end

% Display the percentage of variance for the number of components
pctVar1 = sum(hl(1).ee(1:numComp(1)));
pctVar2 = sum(hl(2).ee(1:numComp(2)));
disp(['% Var, d1 = ' sprintf('%0.1f',pctVar1)]);
disp(['% Var, d2 = ' sprintf('%0.1f',pctVar2)]);

% Now we have two datasets with appropriately* scaled PC scores.

% These just need to be combined...
mzF = repmat(1:numComp(1),[1 numD]);
spF = horzcat(hl.ss);

end


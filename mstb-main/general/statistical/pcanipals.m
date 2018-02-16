function [P,T,varPC] = pcanipals(X,nPCs,tol)
%% pcanipals calculates PC scores and loadings via NIPALS algorithm
%%  Input:   X    - [samples x variables] data matrix, assumed to be mean centered
%            nPCs - number of principal components
%            tol  - tolerance for convergence. Default 1e-6;
%%   Output: T - PC scores [samples x PCs]
%            P - PC loadings [variables x PCs]
%% Author: Kirill Veselkov, Imperial College London 2011.

if nargin < 3
    tol = 10.^-6;
end
    
[nSmpls,nVRbls] = size(X);
T               = NaN(nSmpls,nPCs);
P               = NaN(nVRbls,nPCs);
varX            = nansum(nansum(X.*X));

wb = waitbar(0,'Doing PCA','Name','PCA');

for iPC = 1:nPCs
    [maxnorm,maxnormind] = max(nansum(X.*X));     % find the column with the maximum norm
    t                    = X(:,maxnormind(1)); % select it as the initial score vector
    
    t(isnan(t)) = 0;    
    
    t0                   = t-t;
    
    
    while norm(t-t0) > tol         % iteration to approach the eigenvector direction
        p  = X'*t;      % Project the matrix X onto t in order to find the corresponding loading p
        p  = p/norm(p); % Normalize the loading vector
        t0 = t;         % Save previous t
        t  = X*p;       % Calculate t scores
    end
    X             = X - t*p'; % deflate variance explained by PCs
    T(:,iPC)      = t;        % store score values
    P(:,iPC)      = p;        % store loadings values
    cumvarPC(iPC) = 100*(1-sum(sum(X.*X))./varX); % cumulated variance
    if cumvarPC(iPC)>=100-eps
        fprintf('#%s PCs explain all variance in the data set',num2str(iPC));
        break;
    end
    waitbar(iPC/nPCs,wb,[int2str(iPC) ' / ' int2str(nPCs)]);
end
T = T(:,1:iPC);
P = P(:,1:iPC);
if nPCs >1
    varPC = [cumvarPC(1) diff(cumvarPC(1:iPC))];
else
    varPC = cumvarPC;
end

delete(wb);

return;
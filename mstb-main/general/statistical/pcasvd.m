function [P,T,varPC] =  pcasvd(X,nPCs)
%% pcasvd calculates PC scores and loadings via singular value
%% decomposition
%    Input:    X    - [samples x variables] data matrix, assumed to be mean centered
%%   Author: Kirill Veselkov, Imperial College 2011

if nargin<2
    nPCs = size(X,1); % number of samples
end

[U,D,P] = svd(X./sqrt(size(X,1)-1),'econ');
eigvals = diag(D.^2);
T       = X*P;
rankX   = sum(eigvals>eps); % the number of informative PCs
varPC   = 100*(eigvals./sum(eigvals));

if nPCs<rankX
    rankX = nPCs;
end

T       = T(:,1:rankX); 
P       = P(:,1:rankX); 
varPC   = varPC(1:rankX);
return;
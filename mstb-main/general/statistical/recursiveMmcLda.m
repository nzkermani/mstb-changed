function [P, T, varDV, W] = recursiveMmcLda(X,y,nDVs) 
% LDA: Recursive Linear discriminant analysis 
%  [P, T] = recursiveMmcLda(X, y) 
%  Input: 
%        X - data matrix [observations x variables]
%        y - class labels [observations x 1]
%     nDVs - number of discriminating vectors
%  Output: 
%        P - discriminant feature loadings
%        T - discriminant feature scores  
%    varDV - variance explained by a given factor
%	Reference: 
%       Li, H., Jiang, T. and Zhang, K., Efficient and Robust Feature
%       Extraction by Maximum Margin Criterion. In Proceedings of
%       the Advances in Neural Information Processing Systems 16,
%      (Vancouver, Canada, 2004), MIT Press.
% 
% Author: Veselkov Kirill, Imperial College London  

[nObs,nVar] = size(X);
classLabel  = unique(y); 
nClasses    = length(classLabel); 
meanX       = mean(X);
X           = X - meanX(ones(1,size(X,1)),:); 

if nargin<3
    nDVs = nClasses -1;
end
nDVs  = max(nDVs,2);

% save original data matrix
Xorig = X;
ssqX  = ssq(Xorig);

% Initialize parameters
T     = zeros(nObs,nDVs);
P     = zeros(nVar,nDVs);
varDV = zeros(1,nDVs);

% Need to do dimensionality reduction on the data if the number of
% variables is large
if nVar > (nObs - nClasses) 
    [P1, X] = pcasvd(X,nObs - 1); 
    nVar    = size(X,2);
    pcaDimRed =1;

else
    pcaDimRed = 0;
end

W     = zeros(nVar,nDVs);
Wpca  = 1;

% For each component (DV)
for j = 1:nDVs
    sampleMean = mean(X);
    Sw         = zeros(nVar, nVar);
    Sb         = zeros(nVar, nVar);
    
    % For each of the various classes... This is the one against all
    % principle
    for i = 1:nClasses
        indcs       = (y == classLabel(i));
        nClassSmpls = sum(indcs);
        classMean   = mean(X(indcs,:));
        classMeanSw = classMean(ones(nClassSmpls,1),:);
        classprior  = nClassSmpls/nObs;
        
        % Calculate the within class scatter
        Sw = Sw + (1./nObs).*((X(indcs,:)-classMeanSw)'*(X(indcs,:) - classMeanSw));
        
        % And then the between class scatter
        Sb = Sb + classprior.*(classMean - sampleMean)'*(classMean - sampleMean);
    end
    
    Sw         = (Sw + Sw')/2;
    Sb         = (Sb + Sb')/2;
    
    % Now eig of within - between class scatter
    [w,D] = eigs(Sb-Sw,1,'la',struct('disp',0));
    if w(1)<0; % to be consistent
        w = -w;
    end
    w         = w./norm(w);
    T(:,j)    = X*w;
    W(:,j)    = Wpca*w;
    p         = X'*T(:,j)./(T(:,j)'*T(:,j));
    X         = X - T(:,j)*p'; 
    [wpca, X] = pcasvd(X);
    nVar      = size(X,2);
    Wpca      = Wpca*wpca;
    P(:,j)    = Xorig'*T(:,j)./(T(:,j)'*T(:,j));
    if j>1
        varDV(j) = (1 - ssq(X)./ssqX)*100 - sum(varDV(1:j-1));
    else
        varDV(j) = (1 - ssq(X)./ssqX)*100;
    end
    if abs(D)<eps
        T = T(:,1:nDVs); 
        P = P(:,1:nDVs); 
        varDV = varDV(1:nDVs);
        break;
    end
end

if pcaDimRed == 1  
    W = P1*W;
end
return;
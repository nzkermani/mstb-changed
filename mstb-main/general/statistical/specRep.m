function [ tab ] = specRep(x,y)
%specRep - calculate the spectral reproducibility according to the method
%detailed in the reference. This normalises each spectrum to sum of 1, then
%determines the inner dot product between each spectrum and the group
%average.
%
% INPUT
% x - [m x n], raw spectra of m observations and n variables
% y - {m x 1} or [m x 1], cell array or vector of m rows detailing which
%       unique groups to perform separate calculations
%
% OUTPUT
% tab - table containing the results
%
%
% RUN
% [tab] = specRep(x,y);
%
% Ref: Dill et al., Anal. Bioanal. Chem., 2011, 401, 1949.

% Ensure that matrices are the correct size
if size(y,2) > 1
    y = y';
end
if size(x,1) ~= size(y,1)
    error('Matrices not the same size');
end
    
% Calculate the Euclidean norm of each spectrum
numS = size(x,1);
eucn = zeros(numS,1);
for n = 1:numS
    eucn(n,1) = norm(x(n,:),2);
end

% Normalise each spectrum to a sum of 1
x = bsxfun(@rdivide,x,eucn);

% Determine unique groupings with which to partition the data
[unqG,~,unqI] = unique(y);
numG = numel(unqG);

% Create a results matrix perhaps...
results = zeros(numG,3);

% Loop through each group
for n = 1:numG
    
    % Indices of this group
    ind = unqI == n;
    
    % Determine group average (mean)
    ref = nanmean(x(ind,:),1);
    
    % Determine dot products
    dots = dot(repmat(ref,[sum(ind) 1])',x(ind,:)');
    
    % Now we calculate the mean / std / rsd for these dot products
    results(n,1:2) = [mean(dots) std(dots)];
    results(n, 3 ) = results(n,1) / results(n,2);       
    
end

% Put all in a table
tab = table(unqG,results(:,1),results(:,2),results(:,3),...
    'VariableNames',{'Group','Mean','STD','RSD'});

end


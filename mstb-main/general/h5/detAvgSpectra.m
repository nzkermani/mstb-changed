function [final,file2,group2,idx] = detAvgSpectra(X,file,group,method)
%detAvgSpectra - determine mean/median spectrum for each group in each file

if nargin == 3
    method = 'mean';
end

% Convert file and group to a single class
allClass = classMany2One([file group]);

% Determine unique files
[unqC,~,indC] = unique(allClass);
numC = numel(unqC);

% Empty matrix
final = zeros(numC,size(X,2));
file2 = cell(numC,1);
group2 = cell(numC,1);
idx = zeros(numC,1);

for n = 1:numC
    
    % Indices of this file/group
    fx = indC == n;
    
    % Calculate the average
    switch method
        case 'mean'
            final(n,:) = nanmean(X(fx,:),1);
        case 'median'
            final(n,:) = nanmedian(X(fx,:),1);
    end
    
    % Update the group information
    i = find(fx == 1,1,'first');
    idx(n,1) = i;
    file2(n,1) = file(i,1);
    group2(n,1) = group(i,1);
    
end



end

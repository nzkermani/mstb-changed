function [final,file2,group2,idx] = bootAvgSpectra(X,file,group,method,qty)
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
%final = zeros(numC,size(X,2));
%file2 = cell(numC,1);
%group2 = cell(numC,1);
idx = zeros(numC,1);

% All boot storage?
allB = struct('X',[],'file2',[],'group2',[]);

for n = 1:numC
    
    % Indices of this file/group
    fx = indC == n;
    
    % Calculate the average
    switch method
        case 'mean'
            tmp = bootstrp(qty,@mean,X(fx,:));
            allB(n).X = tmp(randperm(qty,3),:);
        case 'median'
            tmp = bootstrp(qty,@median,X(fx,:));
            allB(n).X = tmp(randperm(qty,3),:);
    end
    
    % Update the group information
    i = find(fx == 1,1,'first');
    idx(n,1) = i;
    
    % What if there is only one sample? Just keep the 1
    if size(allB(n).X,2) == 1
        
        % Then only one sample
        allB(n).X = X(fx,:);
        allB(n).file2 = file(i,1);
        allB(n).group2 = group(i,1);        
    
    else        
        
        % Then there are multiple samples
        allB(n).file2 = repmat(file(i,1),[3 1]);
        allB(n).group2 = repmat(group(i,1), [3 1]);
        
    end
    
end

% Combine all!
final = vertcat(allB.X);
file2 = vertcat(allB.file2);
group2 = vertcat(allB.group2);


end

function [ output_args ] = summaryStat(x,grp)
%summaryStat - determine summary statistics


% Determine groups
[unq,~,ind] = unique(grp);
numG = numel(unq)

res = cell(6,numG+2);
res(1,:) = ['Stat.'; unq; 'Pooled']';

% Mean values...
res(2,1) = {'Mean'};
i = 2;
for n = 1:numG+1
    
    res{i,n+1} = mean(ss(idx == n)

res


end


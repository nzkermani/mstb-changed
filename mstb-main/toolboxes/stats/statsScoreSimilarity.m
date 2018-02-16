function [op,vals] = statsScoreSimilarity(fixed,check)
%statsScoreSimilarity - match the scores in 'check' for consistency against
%those in fixed. Where there are differences, then multiply the appropriate
%axes in check by -1. We need to also have the group indices (numeric) for
%comparison purposes...
%
% James McKenzie, 2017

% Just copy to start with
op = check;

% Determine the number of variables
numV = size(check,2);

% Matrix to store factor values
vals = ones(1,numV);

% Loop through each dimension
for n = 1:numV
    
    % Calculate distance between the points
    par1 = bsxfun(@minus,fixed(:,n),check(:,n));
    
    % Calculate flipped correlation
    par2 = bsxfun(@minus,fixed(:,n),check(:,n) * -1);
    
    % Smallest distance wins...
    if sum(par2) < sum(par1)
        vals(n) = -1;
    end        
    
    
end

% Now just multiply 'check' by 'vals'
op = bsxfun(@times,check,vals);




end


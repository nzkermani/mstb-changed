function [B,W] = oaaMMC(spTrain,groupID,numComps)
% Perform one against all MMC
% 
% INPUTS
%   spTrain     - the [m x n] data matrix
%   groupID     - the [m x 1] vector of class membership (no test set)
%   numComps    - number of DVs to calculate (2 by default)

if nargin == 2
    numComps = 2;
end

% Now we'll begin the OAA approach
[unq,~,ind] = unique(groupID);
numG = numel(unq);

% Somewhere for BETA coefficients
B = zeros(2,numG);       
    
% Matrix for weights
W = zeros(size(spTrain,2),numG);

% Loop
for n = 1:numG
    
    % This is the group that is going up against the others
    tmpIdx = ind == n;
    
    % Run the function...
    [~,T,~,Wi] = recursiveMmcLda(spTrain,tmpIdx,numComps);

    % Single function for the regression coefficients
    [B,W,~] = getRegCoef(n,T,tmpIdx,W,Wi,spTrain,B);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
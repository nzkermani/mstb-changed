function [B,W,T] = getRegCoef(i,T,yi,W,Wi,X,B)
% Get the regression coefficients when using a one against all method

% Correlation
%R = corrcoef([T(:,1),yi]);

% This modification allows us to use the best component rather than using
% only the first component - in PCA, 1st never likely the best.
r = corr(T,yi);
absR = abs(r);
[~,indC] = max(absR);
%indC = 1;


% Flip if negative...based on only the first component (!)
if r(indC) < 0
    W(:,i) = -Wi(:,indC);
else
    W(:,i) = Wi(:,indC);
end

%T = [];

% Scores equal X times weights
T(:,1) = X*W(:,i);

% Get regression coefficients
warning('off','all');
B(:,i) = glmfit(T(:,1),yi, 'binomial', 'link', 'logit');
warning('on', 'all');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

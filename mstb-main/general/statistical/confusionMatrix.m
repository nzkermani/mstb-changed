function [cmat] = confusionMatrix(act,pred,thresh)
% Given assignments of pixels calcualate whether it was done correctly.
% Pred is the Pr matrix from calcScoresSpectra following MMC.  Acutal is a
% list of each pixel's class in a column vector

if nargin == 2
    thresh = 0.80;
end

% Absolute predictions
pred = pred >= thresh;

sa = size(act);
sp = size(pred);

if sa(1) ~= sp(1) 
    error('wrong size');
end

numG = sp(2);

% Extra column is for missed pixels
cmat = zeros(numG,numG+1);

% For each pixel...
for n = 1:sa(1)
    
    % What is the actual class of this pixel?
    clAct = act(n,1);
    
    if clAct ~= 0
        % Was it predicted as anything?
        if sum(pred(n,:)) == 1
            % Single match - what was its class?
            clPred = find(pred(n,:) == 1);
            cmat(clAct,clPred) = cmat(clAct,clPred) + 1;
        else
            cmat(clAct,numG + 1) = cmat(clAct,numG + 1) + 1;
        end
    end    
end

%cmat = 100 * bsxfun(@rdivide,cmat,sum(cmat,2));
cmat(isnan(cmat)) = 0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

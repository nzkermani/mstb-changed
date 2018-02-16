function [T,Pr,Ps] = calcScoresSpectra(allX,meanX,W,B)
% Calculate the scores for spectra, and then their probabilities    

T = bsxfun(@minus,allX,meanX) * W;
nGrps = size(T,2);

% For an actual probability value according to Kirill's code
Pr = zeros(size(T));

for i = 1:nGrps
    %Pr(:,i) = 1 ./ (1 + exp(-(B{i}(1)+B{i}(2)*T(:,i))));
    Pr(:,i) = 1 ./ (1 + exp( -( B(1,i) + B(2,i)*T(:,i) ) ) );
end

% This is to see about non/single/double/multiple classification patterns, 
% and is the same as from SVM so the results are consistent.  A 0.9 
% threshold has been applied throughout 'the toolbox' as the default value.
Ps = Pr >= 0.9;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

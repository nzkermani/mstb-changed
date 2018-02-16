function [ data ] = predictionNormalisation(data)
%predictionNormalisation - check to see what kind of normalisation must be
%performed on the data

numF = size(data,2);
numV = size(data(1).spAl,2);

% Reference spectrum per file
refSpec = NaN(numF,numV);
for n = 1:numF    
    
    % Determine the reference spectrum..
    refSpec(n,:) = nanmean(data(n).spAl(data(n).tobg,:),1);
    
    % Just do TIC normalisation and then quit
    %data(n).spAl = bsxfun(@rdivide,data(n).spAl,nansum(data(n).spAl,2));% * 1e4;
end

mask = sum(refSpec > 0,1) == numF;
for n = 1:numF
    data(n).spAl = data(n).spAl(:,mask);
    data(n).spAl = bsxfun(@rdivide,data(n).spAl,nansum(data(n).spAl,2));
    data(n).mzAl = data(n).mzAl(mask);
end

end


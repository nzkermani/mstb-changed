function [mz,sp] = b2iFilter(dims,bin,thresh,method)
% b2iFilter - James' attempt to make this function more useable to prevent
% iamges of extraordinary size being formed, only to be later shrunk by
% noise filtering functions.  If we filter images here, then we don't need
% to make such a large sp matrix...


switch lower(method)
    case 'withinimage'
        
        % Run function
        [mz,sp] = func1(bin,dims,thresh);
        
    case 'unknown'
        
        % Run other function
        func2(bin,dims);
        
    otherwise
        error('No otherwise');
        
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp] = func1(bin,dims,thresh)
% Main function for binning

% Create a waitbar
wb = waitbar(0,'Filtering');

% How many features are there at the moment?
numF = size(bin,2);
qty = zeros(numF,1);

% The first pass is to identify which of these are going to be useful...
for n = 1:numF
    qty(n,1) = size(bin(n).mz,1);
end

% This is the minimum number of required pixels
minPix = thresh * dims(1) * dims(2);

% These are the bins that we are going to keep
fx = qty >= minPix;

% Create matrix for them
numP = sum(fx);
mz = zeros(numP,1);
sp = zeros(dims(1)*dims(2),numP);

idx = find(fx);

for n = 1:numP
    
    % Index
    i = idx(n);
    
    % Indices
    ii = bin(i).spectra;
    [x,y,~] = ind2sub(dims,ii);
   
    % Convert x|y to single column indices
    i2 = ((y-1)*dims(1)) + x;
    sp(i2,n) = bin(i).intensity;
    
    % Save the mz value
    mz(n,1) = bin(i).centroidmz;
    
    % Update waitbar
    waitbar(n/numP,wb);
   
end

% Reshape into an image
sp = reshape(sp,[dims(1) dims(2) numP]);

% Delete the waitbar
delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function func2(bin,dims)
% Secondary function, not used as far as I am aware

SPnew = zeros(dims);
MZnew = zeros(dims);
for i = 1:nMz
    [temp1,temp2] = ind2sub(dims,bin(i).spectra);
    for l = 1:length(bin(i).spectra)
        SPnew(temp1(l), temp2(l)) = bin(i).centroidspectra;
        MZnew(temp1(l), temp2(l)) = bin(i).centroidmz;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ op ] = cellFreq(ip)
%cellFreq - frequency of elements in a cell, a bit of a text-based
%histogram

% UNique entries
[unq,~,ind] = unique(ip);
numF = numel(unq);
frq = cell(numF,1);

% Loop
for n = 1:numF
    
    fx = ind == n;
    
    frq{n,1} = sum(fx);
    
end

op = [unq frq]

end


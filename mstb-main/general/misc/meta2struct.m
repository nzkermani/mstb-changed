function [new] = meta2struct(meta,heads)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


numF = numel(heads);

for n = 1:numF

    new.(heads{n}) = meta(:,n);
    
end
    

end


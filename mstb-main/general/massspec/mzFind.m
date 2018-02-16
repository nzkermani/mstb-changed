function [ idx ] = mzFind(mzVec,mzQuery,ppmTol)
%mzFind - get indices of m/z values in the main list

if nargin == 2
    ppmTol = 1;
end

numV = numel(mzQuery);

idx = false(size(mzVec));

for n = 1:numV
    
    % Mass tolerance
    df = ppmTol * mzQuery(n) / 1e6;
    %ml = mzFind(n) - df;
    %mh = mzFind(n) + df;
    
    % This will find only a single m/z value if within the range
    [val,tmp] = min(abs(mzVec-mzQuery(n)));    
    if val < df;
        idx(tmp) = true;
    end

end

%disp(num2str(mzVec(idx)))

%disp(['Found ' int2str(sum(idx)) ' / ' int2str(numV) ' peaks']);

end


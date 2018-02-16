function [ fx ] = varsInFiles(sp,files)
%varsInFiles - determine frequency of variables in files given sp matrix
%and list of files

% Determine unique files
[unq,~,ind] = unique(files);
numF = numel(unq);

% Matrix
vec = zeros(numF,size(sp,2));

for n = 1:numF    
    fx = ind == n;    
    vec(n,:) = max(sp(fx,:),[],1) > 0;    
end

% Determine common variables
fx = sum(vec,1) > (numF/2);

end


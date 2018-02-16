function [ new ] = predProbs2Images(pr,cols,thresh)
%predProbs2Images - convert probablity images into something fixed and
%absolute...

% Need prob image (up to 3 layers), appropriate number of colours and a
% threshold

% Check size
sz = size(pr);
if sz > 3
    error('Cannot function');
end

% Reshape and then create new image
pr = reshape(pr,[sz(1)*sz(2) sz(3)]);
new = zeros(sz(1)*sz(2),3);

% Find unambiguously classified pixels
pr = pr >= thresh;
sumC = nansum(pr,2) ~= 1;
pr(sumC,:) = 0;

% Now multiply the columns by 1:sz(3)
cl = bsxfun(@times,pr,[1:sz(3)]);
cl = max(cl,[],2);

% Unique classes and indices
[unq,~,ind] = unique(cl);
if unq(1) == 0
    cols = cat(1,[0 0 0],cols);
end

for n = 1:numel(unq)
    
    % Indices of pixels of this class
    fx = ind == n;
    
    new(fx,:) = repmat(cols(n,:),[sum(fx) 1]);
    
end      
        
new = reshape(new,[sz(1) sz(2) 3]);

end


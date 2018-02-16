function [ output_args ] = detImgOutline(img,flag)
%detImgOutline - calculates the outlines of the discrete values in the
%image

% Only accept single channel images
if size(img,3) > 1
    error('Cannot process multi-channel images');
end

% Fancy plotting it?
if nargin == 2
    flag = true;
else
    flag = false;
end

% Determine unique values
unq = unique(img(:));
if unq(1) == 0
    unq = unq(2:end);
end
numG = numel(unq);

figure; hold on;

% A general place in which to store the results
res = zeros(size(img));

for n = 1:numG
    
    % This is the image mask for this value
    fx = img == unq(n);
    
    % Use the gradient method from...
    [gx,gy] = gradient(double(fx));
    fx((gx.^2 + gy.^2) == 0) = 0;
    
    res = res + (fx * unq(n));
    
    if flag
        [ty,tx] = find(fx == 1);
        scatter(tx,ty,20,'o','filled');
    end
    
end
figure; imagesc(res);

end


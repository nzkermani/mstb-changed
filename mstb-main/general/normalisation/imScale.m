function [ img ] = imScale( img )
%imScale - take an image of n dimensions and scale each to be between 0 and
%1


numF = size(img,3);

for n = 1:numF
    
    tmp = img(:,:,n);
    
    % Make lowest value now equal to 0
    mn = min(tmp(:));
    tmp = tmp - mn;
    
    % Scale top down to 1
    mx = max(tmp(:));
    tmp = tmp ./ mx;
    
    % Replace
    img(:,:,n) = tmp;
    
end
    


end


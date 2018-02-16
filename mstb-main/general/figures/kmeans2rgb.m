function [ new ] = kmeans2rgb(img,cols)
%kmeans2rgb - convert a k-means image of 1---n to an RGB image using the
%cols from a colourmap

% How many unique clusters?
numClust = max(img(:));

% Size of the image
sz = size(img);

% New image
new = zeros([sz 3]);

% If we provide a colourmap, then need to select equally spaced colours
% instead of all of it...
numCol = size(cols,1);
if numCol > numClust    
    cidx = round(linspace(1,numCol,numClust));
    cols = cols(cidx,:);
end

% Now loop through each of the clusters and RGBify them
for n = 1:numClust
    
    fx = img == n;
    
    for r = 1:3
        tmp = new(:,:,r) + (fx * cols(n,r));
        
        
        new(:,:,r) = tmp;
    end

end

end


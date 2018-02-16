function imgErd = myimerode(img,SE)
% myimerode fast image erosion of a binary/grayscale image I 
% with a recangle structure element SE
%  Input:  img - grayscale or binary image [height x widht x 1]
%          SE  - rectangle structural matrix (e.g. [1 1; 1 1])
%  Output: imgErd - eroded image (imgErd)
%% Author: Kirill Veselkov, Imperial College London, 2014

% normalise a structural element (SE) for convolution operation
SE     = SE/sum(SE(:));
% perform image padding for erosion
img    = padImgErs(img,SE);
% convolve the image with the SE
imgErd = conv2(img,SE,'valid');
if length(unique(img))==2
    eps    = 0.001; % define a small constant   
    imgErd = fix(imgErd+eps)==1;
end
end

function img = padImgErs(img,SE)
% matrix padding for errosion
% get the size of the image 
[h,w, ~]  = size(img); 
% compute the size of the stuctural element
[hSE,wSE] = size(SE);
% height and width of the structural element
hSE = hSE - 1;  wSE = wSE - 1;

%% symmetrical matrix padding of the bottom and top parts of the matrix
% (the padding value for the erosion is 1)
if wSE > 1 
    upSE   = floor(wSE/2);  downSE = floor(wSE/2) + rem(wSE,2);
    img    = [ones(h,upSE) img ones(h,downSE)];
    w      = w + upSE + downSE;
elseif wSE == 1
    img = [img ones(h,wSE)]; w = w + wSE;
end

%% symmetrical matrix padding on on the left and right parts of the matrix
if hSE > 1 
    leftSE  = floor(hSE/2);
    rightSE = floor(hSE/2) + rem(hSE,2);
    img     = [ones(leftSE,w); img; ones(rightSE,w)];
elseif hSE == 1
    img     = [img; ones(hSE,w)];
end
end
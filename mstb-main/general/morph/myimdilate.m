function imgDlt = myimdilate(img,SE)
% myimerode fast image dilation of a binary/grayscale image I 
% with a recangle structure element SE
%  Input:  img - grayscale or binary image [height x widht x 1]
%          SE  - rectangle structural matrix (e.g. [1 1; 1 1])
%  Output: imgErd - eroded image (imgErd)
%% Author: Kirill Veselkov, Imperial College London, 2014

% perform image padding for dilation
img    = padImgDlt(img,SE);
% convolve the image with the SE
imgDlt = conv2(img,SE,'valid');
% define a small constant
if length(unique(img))==2
    eps    = 0.001;    
    imgDlt = fix(imgDlt+eps)>=1;
end
end

function img = padImgDlt(img,SE)
% perform
% get the size of the image 
[h,w, ~]  = size(img); 
% compute the size of the stuctural element
[hSE,wSE] = size(SE);
% height and width of the structural element
hSE = hSE - 1;  wSE = wSE - 1;

%% symmetrical zero padding of the bottom and top parts of the matrix
% (the padding value for the dilution is 1)
if wSE > 1 
    downSE = floor(wSE/2);  upSE = floor(wSE/2) + rem(wSE,2);
    img    = [zeros(h,upSE) img zeros(h,downSE)];
    w      = w + upSE + downSE;
elseif wSE == 1
    img    = [zeros(h,wSE) img]; w = w + wSE;
end

%% symmetrical zero padding on on the left and right parts of the matrix
if hSE > 1 
    rightSE = floor(hSE/2);
    leftSE  = floor(hSE/2) + rem(hSE,2);
    img     = [zeros(leftSE,w); img; zeros(rightSE,w)];
elseif hSE == 1
    img     = [zeros(hSE,w); img];
end
end
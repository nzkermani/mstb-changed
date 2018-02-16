function [imRes] = dpnInterpolate(dx)
%dpnInterpolate - take the ion image and interpolate it across the empty
%rows

% First identify the empty / non-empty rows...
fx = sum(dx.img,2) == 0;

% Create the image subset...
imTmp = dx.img(~fx,:);

% Resize the image to full size
imRes = imresize(imTmp,size(dx.img));

end


function [ output_args ] = dpnInterpJSM(img,isInterp,method)
%dpnInterpJSM - method for bilinear and bicubic interpolation

% Provide the [i x j x n] image as fully sized, but with image of isInterp
% values to tell us which values to interpolate.

% Which method to use
switch method
    case 'bilinear'
        neighFunc = @neighBilinear;
        interpFunc = @interpBilinear;
    case 'bicubic'
        neighFunc = @neighBicubic;
        interpFunc = @interpBicubic;
end

% Create a copy of the image...
res = img;

% Size of the image...
sz = size(img)

% Run through each pixel and determine the mask indicating which pixels to
% use for interpolation...
for i = 3:sz(1)-3
    
    % Skip non-interpolated pixels... (actually rows here)
    if ~isInterp(i,1)
        continue;
    end    
    
    % Only run down rows which are to be interpolated
    for j = 3:sz(2)-3
        
        % Find neighbours...
        thisPix = [i j];
        [mask] = neighFunc(thisPix,sz);
        
        % Now interpolate for each ion image...
        for n = 1:sz(3)
            res(i,j,n) = interpFunc(img(:,:,n),mask);            
        end
    end
    
    disp(int2str(i));
    
end
        
        




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask] = neighBilinear(pix,sz)
% Neighbours for bilinear interpolation


% Addresses of the pixels
mask = [pix(1)-1 pix(2)-1; ...
    pix(1)-1 pix(2)+1; ...
    pix(1)+1 pix(1)-1; ...
    pix(1)+1 pix(2)+1];


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask] = neighBicubic(pix,sz)
% Neighbours for bicubic interpolation



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = interpBilinear(img,mask)
% Could just calculate the mean value for this

top = [img(mask(1,1),mask(1,2)) img(mask(2,1),mask(2,2))];
%top = interp1([-1 1],tmp,0);

bot = [img(mask(3,1),mask(3,2)) img(mask(4,1),mask(4,2))];
%bot = interp1([-1 1],tmp,0);

%val = interp1([-1 1],[top bot],0);

val = mean([mean(top) mean(bot)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



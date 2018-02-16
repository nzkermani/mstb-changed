function [ imgn ] = iconProcess( alpha,cdata )
%iconProcess - make an icon good

%close all

%cdata = double(cdata);
%alpha = double(alpha);

% Combine the streams...
img = bsxfun(@plus,cdata,alpha);
%img = cdata;
%figure; imagesc(img);

% Now decide on bg colour...
tr = img(1,end,:); % this may need changing if the image is coloured in the corners...

imgn = NaN(size(img));

% Now run through each pixel, setting to NaN if the same bg
sz = size(img);
bg = zeros(sz(1),sz(2));

for r = 1:sz(1)
    for c = 1:sz(2)
        %img(r,c,:)
        t = img(r,c,:) == tr;
        bg(r,c) = sum(t);
        if bg(r,c) == 3
            img(r,c,:) = 0/0;
        else
            imgn(r,c,:) = 255 - img(r,c,:);
        end
    end
end
imgn = imgn / 255;

return


img2 = bsxfun(@times,imgn,double(cdata)) / 255;

figure; imagesc(img2)


%figure; imagesc(bg);
%figure; imagesc(img2); 
%figure; imagesc(imgn / 255);
% bg = bg == 3;
% nb = ones(sz(1),sz(2));
% nb(bg) = NaN;
% figure; imagesc(nb);
% 
% img2 = img2 / 255;
% img2 = bsxfun(@times,double(img2),nb);








return










% Invert, normalise
%img = 255 - img;
img = bsxfun(@rdivide,double(img),[255 255 255]);

figure; imagesc(img);

% Now find the transparent colour...
tr = img(1,1,:)

% Find those pixels...
bg = bsxfun(@eq,img,tr);
bg = sum(bg,3);
zz = bg == 3;
xx = bg ~= 3;
figure; image(bg)
% Produce NaN mask
bg(zz) = NaN;
%bg(xx) = bg(xx) + 20;

nnn = bg == 1;
img(nnn) = NaN;




figure; image(bg);
figure; image(img);
img2 = bsxfun(@times,double(img),double(bg));
img2 = img2 / 255;
% 
% % Set white to black...
% bg = sum(img2,3);
% ft = ~isnan(bg);




end


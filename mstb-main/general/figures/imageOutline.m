function [ ol ] = imageOutline(img,flag)
%imageOutline - determine the outline of an image, especially for H&E
%images, which will be useful for assessing image coregistration

% Images must be either of three or one dimensions...
sz = size(img,3);
if sz == 3 || sz == 1
else
    img = nansum(img,3);
end

% Need to determine the smoothing factor - depends on optical or MS image
% being used to smooth
if size(img,1) < 200 && size(img,2) < 200
    smoothFactor = 2;
elseif size(img,1) > 500 && size(img,2) > 500
    smoothFactor = 5;
else
    smoothFactor = 3;
end

% Convert to 'gray' by taking the mean
gry = nanmean(img,3);

% Calculate a histogram between integer values (NB requirement to have
% image between 0--256...
if max(gry(:)) <= 1
    gry = gry .* 256;
end
[frq,vals] = hist(gry(:),0:1:256);

% Smooth the histogram...
[sm] = movingWindow(frq,2,'mean');

% Can we find the local minimum somewhere between 200 and the local maximum
% towards the white side of the histogram
mask = vals > 220;
[mx,~] = max(sm(mask));

% Find the maximum...
fx = find(sm == mx);

% Now the minimum between 220 and fx
mask = vals > 200 & vals < fx;
[mx,~] = min(sm(mask));

% This becomes the threshold...
fy = find(sm == mx);
fy = fy(1);

% Determine the outlines of the various bits...
[bgSm,b] = getOL(gry < fy,smoothFactor);

% Determine the size of the outlines - we obvs want the largest one!
szz = cellfun(@max,cellfun(@size,b,'UniformOutput',false));
[~,idx] = max(szz);
ol = b{idx};

% Return from here if no plot requested...
if nargin == 1
    flag = false;
elseif ~islogical(flag)
    flag = false;
end
if ~flag
    return
end

% Now split up the image into 2 parts
imbg = img;
imim = img;
mask = bgSm == 1;
for n = 1:size(img,3)
    
    tmp1 = img(:,:,n);
    tmp2 = img(:,:,n);
    tmp1(mask) = 0;
    tmp2(~mask) = 0;
    
    
    imbg(:,:,n) = tmp1;
    imim(:,:,n) = tmp2;
end
        
figure;

ax(1) = subplot(2,2,1); hold on;
imagesc(img);
plot(ol(:,2),ol(:,1),'LineWidth',2,'Color','r');
xlim([0.5 size(img,2)+0.5]);
ylim([0.5 size(img,1)+0.5]);
set(gca,'YDir','reverse');

ax(2) = subplot(2,2,2);
imagesc(bgSm);

ax(3) = subplot(2,2,3);
imagesc(imbg);

ax(4) = subplot(2,2,4);
imagesc(imim);

linkaxes(ax,'xy');



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bgSm,b] = getOL(img,smFac)
% Take a binary image, smooth it, and then determine the outline(s)

% Now determine the outline of this image. First this will require some
% smoothing, and then the outline can be determined
fh = fspecial('average',smFac);
bgSm = filter2(fh,double(img)) > 0.5;

% Determine boundaries
[b,~,~,~] = bwboundaries(bgSm,8);

% Delete all with a size smaller than 8 connected pixels (from the BG)

bSize = cellfun(@max,cellfun(@size,b,'UniformOutput',false));
chk = bSize >= 10;
b = b(chk);

% Consider deleting the first one too...
try
    if sum(b{1}(1:10,1) == 1) == 10
        b = b(2:end);
    end
catch
    
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

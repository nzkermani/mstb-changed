function [ output_args ] = bg2white(img)
%bg2white - set an image background to white

close all;

% Set final channel to white
bg = img(:,:,3);
bg(bg > 0.5) = NaN;
img(:,:,3) = bg;

% Generic low values to NaN
img(img < 0.1) = NaN;

% Now triple NaN to white
tripNan = sum(isnan(img),3) == 3;
for n = 1:3
    tmp = img(:,:,n);
    tmp(tripNan) = 1;
    img(:,:,n) = tmp;
end

h1 = imagesc(img);
set(h1,'AlphaData',1-tripNan);
set(gca,'XTick',[],'YTick',[],'Box','off','Color','none','Position',[0 0 1 1]);
set(gcf,'Color','none');
return

figure; imagesc(img);
box off
axis off
%set(gca,'Position',[0 0 1 1]);
set(gcf,'color','none');
%set(gca,'color','none');

end


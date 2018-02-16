function [ output_args ] = xxxHistoChange(src,event,fig,man)
%xxxHistoChange - change the histology image transparency

% Early return
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Determine the src, as this function is set off by either the slider or
% the edit box
type = get(src,'Style');
switch type    
    case 'slider'        
        val = get(src,'Value');
        set(man.value,'String',int2str(val));
    case 'edit'        
        val = str2double(get(src,'String'));
        set(man.slide,'Value',val);
end

% Low res optical image?
img = dpn.opt.coreg;
try
    newsz = size(dpn.opt.lowRes);
catch
    newsz = size(dpn.d1.img);
end
img = imresize(img,[newsz(1) newsz(2)]);

% Make a checkerboard image
[sx,sy,~] = size(img);
bs = 5;
p = ceil(sx / bs);
q = ceil(sy / bs);

%al = checkerboard(bs,p,q) * val / 100;
al = ones(sx,sy) * val / 100;

al = linspace(0,val/100,sy);
al = repmat(al,[sx 1]);

figure; hold on;
msimg = dpn.d1.img;
i = imagesc(msimg); colormap(gray);
h = imagesc(img);
set(h,'AlphaData',al);
set(gca,'YDir','reverse');



%set(dpn.fig.ax.ms1(2),'AlphaData',al);





end


function dpnChangeSegImage(src,~,fig,seg,choice,value)
% Change the image on display...

% Get the guidata...
dpn = guidata(fig.fig);

% Determine the appropriate image to be displayed...
nextImg = get(src,'UserData');

% Determine the new data...
new1 = dpn.d1.mva.(choice).img(:,:,nextImg);
if strcmp(dpn.mode,'dual')
    new2 = dpn.d2.mva.(choice).img(:,:,nextImg);
end

% Alter the userdata for the back and forward buttons...
if nextImg-1 > 0
    valBack = nextImg-1;
else
    valBack = size(dpn.d1.mva.(choice).img,3);
end
set(seg.imgPrev,'UserData',valBack);

if nextImg+1 > size(dpn.d1.mva.(choice).img,3)
    valFord = 1;
else
    valFord = nextImg+1;
end
set(seg.imgNext,'UserData',valFord);

% Change the images...
dpnIonImage([],[],dpn.fig.ax.ms1,new1);
if strcmp(dpn.mode,'dual')
    dpnIonImage([],[],dpn.fig.ax.ms2,new2);
end

% Change the text...
set(fig.txtPos,'String',dpn.d1.mva.(choice).label{nextImg});
if strcmp(dpn.mode,'dual')
    set(fig.txtNeg,'String',dpn.d2.mva.(choice).label{nextImg});
end

if size(dpn.d1.mva.(choice).img,3) >= 3 && valFord == 1
    figure; imagesc(imScale(dpn.d1.mva.(choice).img(:,:,[3 2 1])));
    
    if strcmp(dpn.mode,'dual')
        figure; imagesc(imScale(dpn.d2.mva.(choice).img(:,:,[3 2 1])));
    end
end

end
function dpnChangeSegImage2(~,~,fig,seg,choice)
% Change the image on display... this version is going to use the RGB boxes
% for specifying which images should be plotted.  Its objective is similar
% to the original function, but uses the boxes rather than the push
% buttons. However, this version has to be capable of summing multiple
% channels together

% Get the guidata...
dpn = guidata(fig.fig);

% First task: determine which values have been clicked
allVal = get(seg.res,'Value');

% What is the size of the image?
sz = size(dpn.d1.mva.(choice).img);

% Loop through each of the three channels
[img1,img2,img3] = sumImages(dpn,sz,choice,allVal);

% How are the results being displayed?
changeLayoutAndDraw(dpn,img1,img2,img3,choice);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img1,img2,img3] = sumImages(dpn,sz,choice,allVal)

% Create blank images for display purposes
img1 = zeros([sz(1) sz(2) 3]);
if strcmp(dpn.mode,'dual')
    img2 = zeros([sz(1) sz(2) 3]);
    img3 = zeros([sz(1) sz(2) 3]);
else
    img2 = [];
    img3 = [];
end

% Loop through the RED GREEN BLUE
for n = 1:3
    
    % What PCs/LVs/etc were selected by the user?
    switch choice
        case 'kMeans'
            vals = allVal{n};
            vals(vals == 0) = [];
    
        otherwise
            vals = allVal{n} - 1;
            vals(vals == 0) = [];
    end
    
    % Leave blank if nothing selected...
    if isempty(vals)
        continue
    end
    
    % Combine those together...
    img1(:,:,n) = sum(dpn.d1.mva.(choice).img(:,:,vals),3);        
    if strcmp(dpn.mode,'dual')
        img2(:,:,n) = sum(dpn.d2.mva.(choice).img(:,:,vals),3);
        img3(:,:,n) = sum(dpn.fuse.mva.(choice).img(:,:,vals),3);
    end
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeLayoutAndDraw(dpn,img1,img2,img3,choice)
% Determine what needs plotting where

% This function works for both single and dual analysis; thus need to
% ensure that it remains compatible for the two type of analysis
switch dpn.mode
    
    case 'single'
        singleModeLayout(dpn,img1,choice);                

    case 'dual'
        dualModeLayout(dpn,img1,img2,img3,choice);
        
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function singleModeLayout(dpn,img1,choice)
% Function to display the single ion mode image

% Change to dual layout
if strcmp(get(dpn.fig.tb.layout,'State'),'off')
    set(dpn.fig.tb.layout,'State','on');
    desiChangeLayout(dpn.fig.tb.layout,[],dpn.fig);
end

% Plot the classification image
switch choice
    
    case 'kMeans'
        img1 = kmeans2rgb(img1(:,:,1),redbluecmap);
        dpnIonImage([],[],dpn.fig.ax.mv,img1);
        
    otherwise
        dpnIonImage([],[],dpn.fig.ax.mv,img1);
end

% Make a composite to show unclassified pixels
switch choice
    case 'MMC'
        n2 = max(img1,[],3);
        n2 = repmat(n2,[1 1 3]);
        dpnIonImage([],[],dpn.fig.ax.sp,n2);
        set(dpn.fig.ax.sp,'YDir','reverse');

    otherwise
        dpnIonImage([],[],dpn.fig.ax.sp,ones([10 10])*0.5);
end               

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dualModeLayout(dpn,img1,img2,img3,choice)
% Function to display the dual ion ion images. We don't need to change the
% layout for this.

switch choice
    case 'kMeans'
        img1 = kmeans2rgb(img1(:,:,1),redbluecmap);
        img2 = kmeans2rgb(img2(:,:,1),redbluecmap);
        img3 = kmeans2rgb(img3(:,:,1),redbluecmap);
        
    otherwise
        % Do nothing special
end

% Set the ion images
dpnIonImage([],[],dpn.fig.ax.ms1,img1);
dpnIonImage([],[],dpn.fig.ax.ms2,img2);
dpnIonImage([],[],dpn.fig.ax.fu, img3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

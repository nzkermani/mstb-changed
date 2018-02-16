function ihcgui
%ichgui - draw a little gui for processing heatmaps from IHC images

% Draw the various parts...
[fig] = guiDraw;

% Callbacks...
set(fig.load,'Callback',{@imageLoad,fig});
set(fig.doClust,'Callback',{@imageCluster,fig});
set(fig.viewClust,'Callback',{@drawClusterIndividual,fig});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = guiDraw

fig.fig = figure('Name','IHC',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);

% Panel for the control things
fig.pan1 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 0.25 1]);

% Panel for the axes
fig.pan2 = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.25 0 0.75 1]);

% How many axes do we want? 4?
fig.ax(1) = axes('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0.5 0.5 0.5]);
fig.ax(2) = axes('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.5 0.5 0.5 0.5]);
fig.ax(3) = axes('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0 0 0.5 0.5]);
fig.ax(4) = axes('Parent',fig.pan2,...
    'Units','normalized',...
    'Position',[0.5 0 0.5 0.5]);

% Controls...

% Load an image into the workspace
fig.load = uicontrol('Parent',fig.pan1,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.1 0.9 0.8 0.05],...
    'String','Load image');

% Perform cluster analysis
fig.doClust = uicontrol('Parent',fig.pan1,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.1 0.8 0.35 0.05],...
    'String','Cluster');

fig.numClust = uicontrol('Parent',fig.pan1,...
    'Style','edit',...
    'Units','normalized',...
    'Position',[0.55 0.8 0.35 0.05],...
    'String','3');

% Somewhere to display the clusters for selection/visualisation
fig.viewClust = uicontrol('Parent',fig.pan1,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.55 0.6 0.35 0.15],...
    'String','NA');

% Make output...
fig.doOutput = uicontrol('Parent',fig.pan1,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.1 0.1 0.8 0.05],...
    'String','Output');

linkaxes(fig.ax,'xy');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageLoad(~,~,fig)
% Load an image

if ~strcmp(getUser,'jmckenzi')
    imgName = 'IHC-Test1.jpeg';
    imgPath = '/Users/jmckenzi/Dropbox/Imperial/Projects/IHC/';
else    
    defP = '/Users/jmckenzi/Dropbox/Imperial/Projects/IHC/';
    [imgName,imgPath] = uigetfile({'*.png; *.jpeg; *.jpg; *.tiff'},...
        'Select',defP);
end

gd.path = imgPath;
gd.name = imgName;
gd.orig = imread([imgPath imgName]);

% Resize...
gd.orig = imresize(gd.orig,0.25);

% Save as the guidata
guidata(fig.fig,gd);

% Probably best to clear axes in another function and draw this one...
imageDisplay([],[],fig);

% Convert the image
imageConvert([],[],fig);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageDisplay(~,~,fig)

gd = guidata(fig.fig);

for n = 1:numel(fig.ax)
    cla(fig.ax(n),'reset');
end

imagesc(fig.ax(1),gd.orig);
set(gca,'XTick',[],'YTick',[]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageConvert(~,~,fig)
% Convert the image to alternative colour scheme

wb = waitbar(0.5,'Converting');

gc = guidata(fig.fig);

cform = makecform('srgb2lab');
gc.conv = applycform(gc.orig,cform);

% Convert to a double?
gc.conv = double(gc.conv);

guidata(fig.fig,gc);

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageCluster(~,~,fig)

% Get the guidata
gc = guidata(fig.fig);
if isempty(gc)
    return
end

% Convert?
if ~isfield(gc,'conv')
    imageConvert([],[],fig);
    gc = guidata(fig.fig);
end

% Determine the number of clusters requested
k = str2double(get(fig.numClust,'String'));

% Size of the image
sz = size(gc.conv);

% Reshape the image
img = reshape(gc.conv,[sz(1)*sz(2) sz(3)]);

% We only need the 'a' and 'b' parts of the colour scheme
img = img(:,[2 3]);

% k-means clustering
[tmp,~] = kmeans(img,k,...
    'Distance','sqEuclidean',...
    'Replicates',1);

% Set the box to contain the number of clusters...
set(fig.viewClust,'String',int2str((1:k)'),...
    'Value',1);

% Reshape
gc.clust = reshape(tmp,[sz(1) sz(2)]);
gc.numCl = k;

% Save as guidata
guidata(fig.fig,gc);

% Now need to draw the cluster results...
imagesc(fig.ax(2),gc.clust);
set(gca,'XTick',[],'YTick',[]);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawClusterIndividual(~,~,fig)

% Get the guidata
gc = guidata(fig.fig);
if ~isfield(gc,'clust')
    return
end

% Colours
cols = parula(gc.numCl);

% Determine the number of the thing selected in the box
c = get(fig.viewClust,'Value');

mask = gc.clust == c;
tmp = repmat(mask,[1 1 3]);
tmp = bsxfun(@times,double(tmp),reshape(cols(c,:),[1 1 3]));
imagesc(fig.ax(3),tmp);


% What about (tmp) an image from the raw data/converted data?
new = gc.orig;
tmp = new(:,:,1); tmp(~mask) = 0; new(:,:,1) = tmp;
tmp = new(:,:,2); tmp(~mask) = 0; new(:,:,2) = tmp;
tmp = new(:,:,3); tmp(~mask) = 0; new(:,:,3) = tmp;

% Convert the ...

% Perhaps we can consider smoothing the image here...
filt = fspecial('average',10);
tmp = filter2(filt,tmp);

imagesc(fig.ax(4),tmp);
set(gca,'XTick',[],'YTick',[]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

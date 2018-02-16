function [data] = ionScroll_OLD(mz,data)
%histoScroll - given a series of image, scroll through them in a good way
%to enable pseudo 3D visualisation...

% James McKenzie, 2014.

%try    
    % Close existing...
    f0 = findobj('Tag','imageScroll');
    delete(f0);

    % If no argin is provided, then we'll read in from an h5 file...
%     if nargin == 0
%         [data] = readFromH5;
%     end

    % First need to draw the figure...
    [fig] = drawWindow(data,mz);
%return

    % Draw initial
    %sliderReDraw(fig,data);

    % Finally set the callback...
    set(fig.sc,'Callback',{@sliderCallback,fig,data,mz});
    set(fig.fig,'WindowScrollWheelFcn',{@sliderCallback,fig,data,mz});

%catch
    %data = [];
%end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = readFromH5
% Read in the data

[fN,fP] = uigetfile('*.h5','');

% Info of this file...
info = h5info([fP fN]);

% How many datasets / slices within this file?
numD = size(info.Groups,1);

% Empty structure
data = struct('op',[],'mz',[]);

% Now loop...
for n = 1:numD
    data(n).op = h5read([fP fN],['/s' int2str(n) '/op']); 
    tmp        = h5read([fP fN],['/s' int2str(n) '/x' ]);
    data(n).mz = sum(tmp,3);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow(data,mz)
% Draw the gui and its bits

fig.colMap = flipud([64,0,75;...
    118,42,131;...
    153,112,171;...
    194,165,207;...
    231,212,232;...
    247,247,247;...
    217,240,211;...
    166,219,160;...
    90,174,97;...
    27,120,55;...
    0,68,27;...
    0 0 0]);
fig.colMap = fig.colMap ./ 256;

fig.fig = figure('Name','Image Scroll',...
    'Units','normalized',...
    'Position',[0 0 0.5 0.8],...
    'Toolbar','none',...
    'Tag','imageScroll',...
    'Resize','on',...
    'Color',fig.colMap(1,:));


% Draw the axes...
fig.ax1 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 1.0 0.8]);
imagesc(data(:,:,1));
colormap(fig.colMap);
fig.cb1 = colorbar('Location','NorthOutside','FontSize',16,'FontWeight','bold');

% Average spectrum
fig.mz = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.05 0.87 0.9 0.10],...
    'XColor','white');
sz = size(data);
tmp = reshape(data,[sz(1)*sz(2) sz(3)]);
meanSp = nanmean(tmp,1);

hold on;
stem(mz,meanSp,'LineWidth',3,'Color','white');
set(gca,'Color','black',...
    'XColor','white',...
    'YColor','white',...
    'XLim',[mz(1)-1.1 mz(1)+1.1],...
    'YTick',[],...
    'FontSize',16,...
    'FontWeight','bold');
box on;

% Add the red marker to tell us where we are!
fig.mark = scatter(mz(1),meanSp(1),80,'red','o','filled');


% Determine how many images there are...
numImg = size(data,3);
ss = 1 / (numImg - 1);

% Add the clicking scroll bar
fig.sc = uicontrol('Parent',fig.fig,...
    'Style','slider',...
    'Units','normalized',...
    'Position',[0.05 0.94 0.9 0.05],...
    'Min',1,...
    'Max',numImg,...
    'Value',1,...
    'SliderStep',[ss ss]);

% Text box...
fig.tb = uicontrol('Parent',fig.fig,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.33 0.8 0.33 0.05],...
    'String',sprintf('m/z = %0.4f',mz(1)),...
    'BackgroundColor',fig.colMap(1,:),...
    'ForegroundColor','white',...
    'FontWeight','bold',...
    'FontSize',40);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderCallback(~,event,fig,data,mz)
% Update the image in the gui

% Which number to draw now?
imgNo = round(get(fig.sc,'Value'));

% Decide if up/down scroll
if isfield(event,'VerticalScrollCount')
    if event.VerticalScrollCount > 0
        imgNo = min([size(data,3) imgNo + 1]);    
    elseif event.VerticalScrollCount < 0
        imgNo = max([1 imgNo - 1]);
    end
end

% Update the callback functions
set(fig.tb,'String',sprintf('m/z = %0.4f',mz(imgNo)));
set(fig.sc,'Value',imgNo);
set(fig.mz,'XLim', [mz(imgNo)-1.1 mz(imgNo)+1.1]);

% Update the red marker
meanVal = data(:,:,imgNo);
meanVal = nanmean(meanVal(:));
set(fig.mark,'XData',mz(imgNo),'YData',meanVal);

% Call the draw function...
sliderReDraw(fig,data);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderReDraw(fig,data)

% Update text box
imgNo = get(fig.sc,'Value');

ih = get(fig.ax1,'Children');
set(ih,'CData',data(:,:,imgNo));

% Draw the image...
%axes(fig.ax1);
%imagesc(data(:,:,imgNo));
%set(gca,'XTick',[],'XTickLabel',[],...
%    'YTick',[],'YTickLabel',[]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

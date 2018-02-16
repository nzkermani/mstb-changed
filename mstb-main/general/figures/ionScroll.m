function ionScroll(mz,data)
%histoScroll - given a series of image, scroll through them in a good way
%to enable pseudo 3D visualisation...

% Sort the data according to m/z
[~,idx] = sort(mz,'ascend');

mz = mz(idx);
data = data(:,:,idx);

% Determine the Matlab version
vrs = ver('Matlab');
if strcmp(vrs.Release,'(R2014a)')
    vrs = 2014;
elseif strcmp(vrs.Release,'(R2016a)')
    vrs = 2016;
end
    

% James McKenzie, 2014.

%try    
    % Close existing...
    f0 = findobj('Tag','histoScroll');
    delete(f0);

    % If no argin is provided, then we'll read in from an h5 file...
%     if nargin == 0
%         [data] = readFromH5;
%     end

    % First need to draw the figure...
    [fig] = drawWindow(data,mz);

    % Draw initial
    sliderReDraw(fig,data);

    % Finally set the callback...
    set(fig.sc,'Callback',{@sliderCallback,fig,data,mz,vrs});
    set(fig.fig,'WindowScrollWheelFcn',{@sliderCallback,fig,data,mz,vrs});

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

fig.fig = figure('Name','Histo Scroll',...
    'Units','normalized',...
    'Position',[0.1 0.1 0.8 0.8],...
    'Toolbar','none',...
    'Tag','histoScroll',...
    'Resize','on');

hold on;

% Draw the axes...
fig.ax1 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0.1 1.0 0.90]);

% Determine how many images there are...
numImg = size(data,3);
ss = 1 / (numImg - 1);

% Add the clicking scroll bar
fig.sc = uicontrol('Parent',fig.fig,...
    'Style','slider',...
    'Units','normalized',...
    'Position',[0.1 0 0.8 0.05],...
    'Min',1,...
    'Max',numImg,...
    'Value',1,...
    'SliderStep',[ss ss]);

% Text box...
fig.tb = uicontrol('Parent',fig.fig,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.4 0.05 0.2 0.02],...
    'String',sprintf('m/z = %0.4f',mz(1)));

fig.tb2 = uicontrol('Parent',fig.fig,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.8 0.05 0.2 0.02],...
    'String','1');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderCallback(~,event,fig,data,mz,vrs)
% Update the image in the gui

% Which number to draw now?
imgNo = round(get(fig.sc,'Value'));

% Various things for different versions...
if vrs <= 2014
    if isfield(event,'VerticalScrollCount')
        if event.VerticalScrollCount > 0
            imgNo = min([size(data,3) imgNo + 1]);    
        elseif event.VerticalScrollCount < 0
            imgNo = max([1 imgNo - 1]);
        end
    elseif isfield(event.Source,'Value')
        imgNo = event.Source.Value;
    end
    
elseif vrs >= 2015 % not sure if 2015 will work with this
    
    switch event.EventName
        
        case 'WindowScrollWheel'            
            if event.VerticalScrollCount > 0
                imgNo = min([size(data,3) imgNo + 1]);
            elseif event.VerticalScrollCount < 0
                imgNo = max([1 imgNo - 1]);
            end
            
        case 'Action'            
            imgNo = round(event.Source.Value);
            
        otherwise
            disp('There is no otherwise');
    end
    

end

    



% % Decide if up/down scroll
% if isfield(event,'VerticalScrollCount')
%     if event.VerticalScrollCount > 0
%         imgNo = min([size(data,3) imgNo + 1]);    
%     elseif event.VerticalScrollCount < 0
%         imgNo = max([1 imgNo - 1]);
%     end
% elseif ~isempty(event.VerticalScrollCount)
%     % This is the adjustment required for version 2016a
%     if event.VerticalScrollCount > 0
%         imgNo = min([size(data,3) imgNo + 1]);    
%     elseif event.VerticalScrollCount < 0
%         imgNo = max([1 imgNo - 1]);
%     end
% elseif isfield(event.Source,'Value')
%     imgNo = event.Source.Value;
% end

% Update the callback functions
set(fig.tb,'String',sprintf('m/z = %0.4f',mz(imgNo)));
set(fig.sc,'Value',imgNo);
set(fig.tb2,'String',imgNo);

% Call the draw function...
sliderReDraw(fig,data);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderReDraw(fig,data)

% Update text box
imgNo = get(fig.sc,'Value');

% Draw the image...
axes(fig.ax1);
imagesc(data(:,:,imgNo));
set(gca,'XTick',[],'XTickLabel',[],...
    'YTick',[],'YTickLabel',[]);

colormap(parula);

% axes(fig.ax2);
% imagesc(data(imgNo).mz); colormap(gray);
% set(gca,'XTick',[],'XTickLabel',[],...
%     'YTick',[],'YTickLabel',[]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

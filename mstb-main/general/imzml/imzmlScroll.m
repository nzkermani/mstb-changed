function [ data ] = imzmlScroll( ip )
%imzmlScroll - read in an imzml file and scroll through the images with a
%defined ppm/Da window

% Define a few constants
res = 0.001;

% What to do with the imports?
doImport = false;
if nargin == 0
    
    % Ask user for the file
    defP = '/Volumes/JSM/DB/';
    if ~exist(defP,'dir')
        defP = pwd;
    end
    
    % Get the imzML file
    [a,b,~] = uigetfile({'*.imzML'},'Select File',defP);
    file = [b a];    
    doImport = true;
    
elseif nargin == 1
    
    % If it is text, check file exists
    if ischar(ip)
        if exist(ip,'file')
            doImport = true;
        else
            error('File does not exist');
        end        
    elseif ~iscell(ip)
        error('Incorrect input format');        
    end
end

% Here is where we could import the imzML file into a cell format
if doImport    
    data = imzmlRawExtract(file);
else
    data = ip;
    clear ip;
end

% Can we determine the mz vector
[mzRange] = detmzvector(data);

% Now we need to draw the window, with two axes and some sliders. Should be
% able to mostly copy ionScroll function
[fig] = drawWindow(mzRange,res,size(data));

% Set callback
set(fig.sc,'Callback',{@sliderCallback,fig,data,mzRange,res});
set(fig.fig,'WindowScrollWheelFcn',{@sliderCallback,fig,data,mzRange,res});



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz] = detmzvector(data)
% Get the mz vector from small to large

lo = NaN(size(data));
hi = NaN(size(data));

for p = 1:size(data,1)
    for q = 1:size(data,2)
        lo(p,q) = min(data{p,q}(:,1));
        hi(p,q) = max(data{p,q}(:,1));
    end
end

mz = [floor(min(lo(:))) ceil(max(hi(:)))];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawWindow(mzRange,res,sz)
% Draw the gui and its bits

fig.fig = figure('Name','Histo Scroll',...
    'Units','normalized',...
    'Position',[0.1 0.1 0.8 0.8],...
    'Tag','imzmlRawExtract',...
    'Resize','on');


% Draw the axes...
fig.ax1(1) = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.05 0.2 0.4 0.75]);

fig.ax1(2) = imagesc(rand(sz));
fig.ax1(3) = colorbar;

fig.ax2(1) = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.55 0.2 0.4 0.75]);

fig.ax2(2) = imagesc(rand(sz));
fig.ax2(3) = colorbar;

% Slider steps
df = mzRange(2) - mzRange(1);
largeStep = 1 / (df - 1);
smallStep = 1 / ((df / res) - 1);

% Add the clicking scroll bar
fig.sc = uicontrol('Parent',fig.fig,...
    'Style','slider',...
    'Units','normalized',...
    'Position',[0.05 0 0.6 0.1],...
    'Min',1,...
    'Max',df,...
    'Value',1,...
    'SliderStep',[smallStep largeStep]);

% Text box...
fig.tb = uicontrol('Parent',fig.fig,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.4 0.05 0.2 0.02],...
    'String',sprintf('%0.4f',mzRange(1)),...
    'FontSize',16);

% Text edit box for ppm tolerance
fig.ppm = uicontrol('Parent',fig.fig,...
    'Style','edit',...
    'Units','normalized',...
    'Position',[0.75 0.05 0.05 0.05],...
    'String','5',...
    'FontSize',16);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sliderCallback(~,event,fig,data,mz,res)
% Update the image in the gui

% Which number to draw now?
imgNo = get(fig.sc,'Value');

% Various things for different versions...  
switch event.EventName
    

    case 'WindowScrollWheel'            
        if event.VerticalScrollCount > 0
            imgNo = min([mz(2)-mz(1) imgNo + (100 * res)]);
        elseif event.VerticalScrollCount < 0
            imgNo = max([1 imgNo - (100 * res)]);
        end
        
    case 'Action'            
        imgNo = event.Source.Value;

    otherwise
        disp('There is no otherwise');
end

% CONVERT imgNo into an m/z value
fx = (imgNo-1) + mz(1);
fx = round(fx / res) * res;

% Determine the ppm value
ppm = str2double(fig.ppm.String);
tol = ppm * fx / 1e6;

% Print to the box
str = [sprintf('%0.4f',fx-tol) ' --> ' sprintf('%0.4f',fx) ' <-- ' sprintf('%0.4f',fx+tol)];
set(fig.tb,'String',str);

% Set the value (esp. important for scroll wheel)
set(fig.sc,'Value',imgNo);

% Call the draw function...
imageDraw(fig,data,fx,tol);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageDraw(fig,data,mz,tol)

% Now need to extract the image
spImg = NaN(size(data));
mzImg = NaN(size(data));
for p = 1:size(data,1)
    for q = 1:size(data,2)-1
        
        
        mask = data{p,q}(:,1) >= (mz-tol) & data{p,q}(:,1) <= (mz+tol);
        
        tmp = data{p,q}(mask,:);
        
        
        if ~isempty(tmp)        
            
            [~,b] = max(tmp(:,2));
            mzImg(p,q) = tmp(b,1);
            spImg(p,q) = tmp(b,2);
            
        end
        
    end
end

% Can we determine the median m/z value for this image?
medMZ = nanmedian(mzImg(:));

% Convert mzImg to ppmImg
mzImg = 1e6 * (mzImg - mz) / mz;

set(fig.ax1(2),'CData',spImg);
set(fig.ax2(2),'CData',mzImg);

title(fig.ax2(1),sprintf('%0.4f',medMZ));%,'FontSize',16);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

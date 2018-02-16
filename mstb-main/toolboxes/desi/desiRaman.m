function desiRaman(src,event,fig)
%desiRaman - upload a Raman sample from Mads/Andrea into the toolbox and
%use the image masks as well.
%
% It requires the Raman images to be loaded into the matlab workspace.

% Reset the toolbox to make it ready to recieve new files

% Close annotation / manipulation windows
f0 = findobj('Tag','manipulation');
close(f0);
f0 = findobj('Name','annotation');
close(f0);

% Change the layout to just the two window view...
xxxEnforceLayoutChange(fig,'single','off');

% Clear the axes... but what if we have annotation boxes?
sz = size(get(fig.ax.opt(2),'CData'));
set(fig.ax.opt(2),'CData',ones(sz));
sz = size(get(fig.ax.ms1(2),'CData'));
set(fig.ax.ms1(2),'CData',ones(sz));
f0 = findobj('Type','patch');
delete(f0);
f0 = findobj('Type','scatter');
delete(f0);

% Need to turn the grid off too...
set(fig.tb.grid,'State','off');
desiGrid(fig.tb.grid,[],fig)

% Load the Raman variables
[x,optImg,flag] = getRamanData;
if ~flag
    error('Some issue with the Raman data');
end

% Now here we should consider performing norm2 on the data
n2n = zeros(size(x,1),size(x,2));
for n = 1:size(x,1)
    for r = 1:size(x,2)
    n2n(n,r) = norm(squeeze(x(n,r,:)),2);
    end
end
x = bsxfun(@rdivide,x,n2n);
    

% Now save the results into the structure as required.
file.dir = '/Users/jmckenzi/Desktop/';
file.nam = 'Unknown-Raman';
file.ext = 'mat';
dpn.file = file;

opts.raman = [];
dpn.opts = opts;

% Try to get the wavenumbers from the workspace if they have been loaded.
% Note that we only have the background subtracted ones.
if size(x,3) == 496
    try
        dpn.d1.mz = evalin('base','axis');
    catch
        dpn.d1.mz = 1:size(x,3);
    end
else          
    dpn.d1.mz = 1:size(x,3);
end

% This is the spectral matrix - it may require unit vec normalisation to
% make it a lot nicer to work with
dpn.d1.sp = x;

dpn.fig = fig;
dpn.defP = ['/Users/jmckenzi/Dropbox/Imperial/Publications/'...
    'Myelin Brain DESI Raman/Raman Data/'];

dpn.opt.coreg = optImg;
dpn.opt.orig  = optImg;

dpn.mode = 'single';

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

% Add the ion image
% Draw the optical image...
set(fig.ax.opt(2),'CData',optImg);
set(fig.ax.opt(1),...
    'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(optImg,2)],...
    'YLim',[1 size(optImg,1)],...
    'TickLength',[0 0]);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,optImg,flag] = getRamanData

flag = true;

% Determine the list of variables in base workspace
vars = evalin('base','who');

% Create a list box to define which 
[a,b] = listdlg('PromptString','Select a Raman Variable',...
    'SelectionMode','single',...
    'ListString',vars);

if b == 0
    flag = false;
    return
end

% Get the variable, and determine size
x = evalin('base',vars{a});
sz = size(x);

% Now we do the same produre, but ask for multiple files to be the mask for
% the image.
[a,b] = listdlg('PromptString','Select a Raman Mask',...
    'SelectionMode','multiple',...
    'ListString',vars);

if b == 0
    flag = false;
    return
end

% Determine all sizes
numMask = numel(a);
masks = cell(numMask,1);
maskSz = zeros(numMask,2);
for n = 1:numMask
    masks{n,1} = evalin('base',vars{a(n)});
    maskSz(n,:) = [size(masks{n,1},1) size(masks{n,1},2)];
end

% Check that there is perfect agreement between sizes of image and mask
chkFun = bsxfun(@eq,maskSz,sz(1:2));
if sum(chkFun(:)) ~= numel(chkFun)
    flag = false;
    return
end

% So now create a summed ion image of the various masks
optImg = zeros(maskSz(1,:));
for n = 1:numMask    
    optImg = optImg + masks{n} * n;
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
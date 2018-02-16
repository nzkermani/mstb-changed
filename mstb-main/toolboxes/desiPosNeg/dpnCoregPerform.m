function dpnCoregPerform(~,~,fig,cor)
%dpnCoregPerform - actually run the coregistration function to rende the
%images together

% Get the guidata
dpn = guidata(fig.fig);

% Determine which ion mode is being coregistered
tmp = get(cor.listMS,'String');
mode = get(cor.listMS,'Value');
mode = tmp{mode};
switch mode
    case 'Positive'
        dx = dpn.d1;
        ax = [dpn.fig.ax.ms1 dpn.fig.ax.ms2];
    case 'Negative'
        dx = dpn.d2;
        ax = [dpn.fig.ax.ms2 dpn.fig.ax.ms1];
    case 'Otsu'
        dx = dpn.d1;
        ax = [dpn.fig.ax.opt dpn.fig.ax.ms1];
end

% Determine if we need to do 1 iteration
skipCoreg = get(cor.skip,'Value');
if skipCoreg
    nIter = 1;
else
    nIter = 1000;
end

% May want to remove the annotations by force as well...
if isfield(dpn,'anno')
    
    f0 = findobj('Type','patch');
    delete(f0);
        
    dpn = rmfield(dpn,'anno');
end

% Find the two pixel images for the coregistration...
try
    alOpt = double(dpn.opt.tobg);
    alMSI = double(dx.tobg);
    %reset(ax(2));
catch
    disp('Cannot find the current pixels to be aligned');
    return
end

% Delete the green box...
try
    delete(dpn.alRes.haffine);
catch
end

% Initialise the coregistration function by drawing stuff in the window
dpn.alRes = setOP2MSinBW(alOpt,[],zeros(size(alOpt)),ax);

% Now we need to start the actual co-registration process...
stWarp  = [0 0 0; 0 0 0];

colormap gray

% Perform the affine transformation
[dpn.alRes.rms_error,...
    dpn.alRes.rec_pnts,...
    dpn.alRes.error_img,...
    dpn.alRes.BWaligned,...
    dpn.alRes.usedIters] = affine_ia2(alOpt,...
                                alMSI,stWarp,nIter,...
                                [dpn.alRes.haffine ax(2)]);

% Delete the green box...
delete(dpn.alRes.haffine);
                            
% Generate the new optical image based on the affine transformation
[dpn.opt.coreg] = genNewOptImg(dpn,dpn.opt.orig);
                            
% Save this information to the structure...
guidata(fig.fig,dpn);

% Update the plots to show all the right stuff again...
set(fig.ax.opt(2),'CData',dpn.opt.coreg,...
    'XData',[1 size(dpn.opt.coreg,2)],...
    'YData',[1 size(dpn.opt.coreg,1)]);
set(fig.ax.opt(1),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.coreg,2)],...
    'YLim',[1 size(dpn.opt.coreg,1)]);

% Calculate the grid (again) to adapt to the (new) co-registration
desiDetermineGrid([],[],fig);

% Update both images...
dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);
if strcmp(dpn.mode,'dual')
    dpnIonImage([],[],fig.ax.ms2,dpn.d2.img);
end

% Change the text labels in the boxes
set(cor.listOP,'Value',1);
set(cor.listMS,'Value',1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alRes = setOP2MSinBW(opobject,~,errobject,ax)

warp_p = [0 0 0; 0 0 0];

% Template size
[h,w] = size(opobject);

% Template verticies, rectangular 
% [minX minY; minX maxY; maxX maxY; maxX minY]
tmplt_pts = [1 1; 1 h; w h; w 1]';

% Initialise the figure
colormap(gray);

% This is for the green box to be drawn on...
axes(ax(3));
hold on;

M           = [warp_p; 0 0 1];
M(1,1)      = M(1,1) + 1;
M(2,2)      = M(2,2) + 1;
warp_pts    = M * [tmplt_pts; ones(1, size(tmplt_pts,2))];

% Draw the green box...
alRes.haffine = plot([warp_pts(1,:) warp_pts(1,1)], ...
    [warp_pts(2,:) warp_pts(2,1)],...
    'g-',...
    'LineWidth',3); 


% Error image, this goes in the final axes
set(ax(2),'CData',errobject);
set(ax(1),'XTick',[],'YTick',[],...
    'XLim',[1 size(errobject,2)],...
    'YLim',[1 size(errobject,1)]);
    
%axes(ax(1));
%alRes.herror = imagesc(errobject); 

alRes.warp_pts = warp_pts;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newOpt] = genNewOptImg(dpn,img)
% Generate the optical image using the transformation matrix

% Redefine the x/y grid?
[rowA,colA,~]   = size(img);
[rowB,colB,~]   = size(dpn.alRes.BWaligned);
xgrid           = rowA/rowB : rowA/rowB : rowA;
ygrid           = colA/colB : colA/colB : colA;
[xi,yi]         = meshgrid(xgrid,ygrid);

% Rec points, now having been corrected for the grid sizes
rec_pnts(1,:)   = dpn.alRes.rec_pnts(1,:) * yi(1); 
rec_pnts(2,:)   = dpn.alRes.rec_pnts(2,:) * xi(1);

% The original image to be redrawn...
[h,w,~]         = size(img);

% Template points
tmplt_pts       = [1 1 1; 1 h 1; w h 1; w 1 1]';

% The affinity warp matrix
affityWarp      = rec_pnts * pinv(tmplt_pts);

% Transform...
newOpt = NaN(size(img));
for r = 1:3
    newOpt(:,:,r) = warp_b(double(img(:,:,r)),affityWarp);
end
newOpt = newOpt ./ 256;

% Now try to be a little smarter with the background - we don't want it to
% be black, but closer to the white of the background. Use the first row to
% determine the generic background colour...
tmp1 = newOpt(1,:,:);
val1 = median(tmp1,2);
tmp2 = newOpt(:,1,:);
val2 = median(tmp2,1);
bgCol = max([val1(:) val2(:)],[],2);

mask = sum(newOpt,3) == 0;
for n = 1:3
    mask2 = mask * bgCol(n);       
    newOpt(:,:,n) = bsxfun(@plus,newOpt(:,:,n),mask2);    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wimg = warp_b(img, M, dst)
% Warp the image img according to the provided rectangle

if nargin < 3
    [h,w] = size(img);
    % Template verticies, rectangular 
    % [minX minY; minX maxY; maxX maxY; maxX minY]
    dst = [1 1; 1 h; w h; w 1]';
end

wimg = quadtobox(img, dst, M, 'bilinear');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

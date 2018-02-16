function [mz,x,img,opt,allAnno,imz] = desiH5LoadFunction(file)
%desiH5Load - load an existing h5 file into the desi toolbox

% Define the name from the 'file' input
fn = [file.dir file.nam];

% Create a waitbar
%wb = waitbar(0,'Loading an H5 file');

% Please try to read the imzML filename from which the raw data was
% acquired. This is stored (usually) in the h5 file, as "paths"
imz.dir = h5readatt(fn,'/paths','path');
imz.nam = h5readatt(fn,'/paths','name');

% Read in the main parts of the file (mz,x,images)
[mz,x,img,opt] = readMainParts(fn);

% Get the annotation parts
[names,annoOpt,annoMS] = readAnnotations(fn);
annoMS = round(annoMS);

% Really want to convert the annotations on the optical image to match the
% grid. The grid needs to be determined which is going to make things quite
% hard...
[tmp.fig.grid.opx,tmp.fig.grid.opy,...
    tmp.fig.grid.msx,tmp.fig.grid.msy] = xxxGridCalculate(opt.coreg,img);

newOpt = zeros(size(annoOpt));
for n = 1:size(annoOpt,1)
    
    c1 = annoMS(n,[1 2]) + [1 2];
    c2 = annoMS(n,[3 4]) + [1 2];
    
    if c1(2) > size(img,2)
        c1(2) = size(img,2);
    end
    if c2(2) > size(img,1)
        c2(2) = size(img,1);
    end
    
    newOpt(n,[1 2]) = tmp.fig.grid.opx(c1);
    newOpt(n,[3 4]) = tmp.fig.grid.opy(c2);
    
end
annoOpt = newOpt;

% Match up the colours...
[newCols,newNames,usd] = colourMatching(fn,names);

% Make the annotation table
[allAnno] = makeAnnotationTable(annoOpt,annoMS,newCols,newNames,usd);

return

% Create a guidata structure, and then save the results, and then display
% them in the right boxes...
dpn.file = file;

% This bit is a little unknown at the moment. Presumably they need to be
% specific to this toolbox...
dpn.opts = [];

% Copy the mz and sp vector/matrices
dpn.d1.mz = mz;
dpn.d1.sp = x;
dpn.d1.img = img;

% Coregistered images...
dpn.opt.coreg = opt.coreg;

% These were provided
dpn.fig = fig;
dpn.defP = defP;

% And the annotations
dpn.anno = allAnno;

% Mandatory for single non PS data
dpn.mode = 'single';

% From here we could do the dpnMatlabLoad function, just skipping the
% reading in functionality part...
dpnMatlabLoad(dpn,dpn.fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,x,img,opt] = readMainParts(fn)

% Load the mz vector
mz = h5read(fn,'/mz');

% Load the sp matrix
x  = h5read(fn,'/X');

% Check for log transformed data
minX = min(x(:));
maxX = max(x(:));
if minX > 0 || maxX < 100    
    x = exp(x) - exp(minX);    
end


img = h5read(fn,'/Xsum');

% Load the optical image (orig,lowRes,gray,tobg,coreg)
opt.orig    = h5read(fn,'/opImage');
opt.lowRes  = h5read(fn,'/lropia');
opt.gray    = [];
opt.tobg    = h5read(fn,'/Xsum');
opt.coreg   = h5read(fn,'/opimage_aligned');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [names,annoOpt,annoMS] = readAnnotations(fn)

% Read the annotation names...
tmp = h5info(fn,'/tissue_id');
numA = size(tmp.Attributes,1) - 1;
names = cell(numA,1);
for n = 1:numA
    names{n,1} = h5readatt(fn,'/tissue_id',int2str(n));
end

% If we have more than 15 annotations, then the toolbox cannot handle the
% file.  Would need to support more colours
if numA > 15
    error('Cannot load file with more the 15 annotation classes');
end


% Get the annotation locations
tmpX = h5read(fn,'/opImage_xrect');
tmpY = h5read(fn,'/opImage_yrect');
annoOpt = [tmpX(:,1:2) tmpY(:,[1 3])];

tmpX = h5read(fn,'/msImage_xrect');
tmpY = h5read(fn,'/msImage_yrect');
annoMS = [tmpX(:,1:2) tmpY(:,[1 3])];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newCols,newNames,usd] = colourMatching(fn,names)
% Match up the colours to those in this toolbox.  Closest matches.

% Colours of the groups defined above - this is required for matching the
% things together...
groupCols = h5read(fn,'/groupColors');

% Colours of each annotation
cols = h5read(fn,'/userROIcolours');

numA = size(groupCols,1);

% We need to match up the colours from the H5 file to those 15 in the
% toolbox. Look at similarity 
tbCols = [215 25 28; 253 174 97; 255 255 191; 166 217 106; 26 150 65;...
    166 97 26; 223 194 125; 245 245 245; 128 205 193; 1 133 113;...
    208 28 139; 241 182 218; 146 197 222; 5 113 176; 64 64 64] / 256;

colMatch = zeros(size(tbCols,1),numA);
for n = 1:numA
    
    % Colour similarity...
    csim = bsxfun(@minus,tbCols,groupCols(n,:));
    
    % Find the closest match
    colMatch(:,n) = sum(abs(csim),2);
end

% Determine the lowest values in each column
lowVal = min(colMatch,[],1);
[~,idx] = sort(lowVal);

% Vector to ensure that the colours are not used multiple times
usedCols = zeros(size(tbCols,1),1);

% Find the best colours...
newGroupCols = zeros(size(groupCols));
newCols = zeros(size(cols));
usd = zeros(size(groupCols,1),1);
newNames = cell(size(groupCols,1),1);
for n = 1:numA
    
    % Get the column and set values to 1 if the colour has been taken
    tmp = colMatch(:,idx(n));
    mask = usedCols == 1;
    tmp(mask) = 10;
    
    [~,cidx] = min(tmp);
    usedCols(cidx,1) = 1;
    newGroupCols(idx(n),:) = tbCols(cidx,:);
    
    % Replace the values in 'cols' with these values
    fx = bsxfun(@eq,cols,groupCols(idx(n),:));
    fx = sum(fx,2) == 3;
    newCols(fx,1) = newGroupCols(idx(n),1);
    newCols(fx,2) = newGroupCols(idx(n),2);
    newCols(fx,3) = newGroupCols(idx(n),3);    
    
    usd(fx,1) = cidx;
    
    newNames(fx,1) = names(idx(n));

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allAnno] = makeAnnotationTable(annoOpt,annoMS,newCols,newNames,usd)

% Anonymous function for coloured cell
colorgen = @(color,text) ['<html><table border=0 width=400 bgcolor=',...
    color,'><TR><TD>',text,'</TD></TR> </table></html>'];

% Put the annotation information in the space that we expect to find it in
% this toolbox
numAnno = size(annoOpt,1);
allAnno = cell(numAnno,9);
for n = 1:numAnno
    
    % Generate the new colour
    col2 = colorgen(rgb2hex(newCols(n,:)),int2str(usd(n,1)));
    
    % Put the stuff in the right places...
    allAnno(n,:) = {[0 0 0] usd(n) newCols(n,:) col2 newNames{n,:} annoOpt(n,1:2) annoOpt(n,3:4) annoMS(n,1:2) annoMS(n,3:4)};
    
    % Place the important information into a cell
    %tmpAnno = {[0 0 0] usd col col2 ['ID ' int2str(usd)] xv yv msX msY};

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
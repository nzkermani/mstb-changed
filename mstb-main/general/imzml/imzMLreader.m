function [ dat,op ] = imzMLreader(varargin)
% imzMLreader - read in the data from an imzML file.
%
% How to use this marvellous function?
% Type the following to make it work, you'll need to change the things that
% appear between <> 
%
% imzMLreader('folder',<a>,'file',<b>,'mz',<c>,'regions',<d>,'conc',<e>);
%
% <a>, folder where the file is
% <b>, name of the file (imzML)
% <c>, m/z range to be summed, e.g. [260.16 260.17]
% <d>, number of regions to draw, e.g. 5
% <e>, optional, the concs of the regions, e.g. [1 5 10 50 100]


% Get the addition stuff
[opts] = readArgsData(varargin);

% What is the file extension?
ext = fileExtension(opts.fN);

% Get the stuff from the file
switch lower(ext)
    case 'imzml'
        [op] = getMZvectors(opts.fP,opts.fN,opts);
    case 'raw'
        [op] = getFromRawFile(opts.fP,opts.fN,opts);
end

% Size of the grid
sz = size(op.grid);
roi = zeros(sz(1),sz(2),opts.numR);
dat = struct('reg',[],'val',[],'sum',[]);

% Draw the image and allow roipoly for drawing the thing
figure('Units','normalized','Position',[0.25 0.25 0.5 0.5]); 
fig = imagesc(log(op.grid));
title([sprintf('%0.4f',opts.mz(1)) ' - ' sprintf('%0.4f',opts.mz(2)) ' Th'],...
    'FontSize',22,'FontWeight','bold');
clab = colorbar;
ylabel(clab,'ln(intensity)','FontSize',16,'FontWeight','bold');
axis image

for n = 1:opts.numR

    % Draw the regions
    roi(:,:,n) = roipoly;
        
    % Gather data points from the region
    idx = roi(:,:,n) == 1;
    dat(n).val = op.grid(idx);
    dat(n).sum = sum(dat(n).val);
    dat(n).reg = repmat(opts.conc(n),sum(idx(:)),1);
end

% Output to a single text file...
dat2txt(opts.fP,opts.fN,dat,opts);

% Concat for box plot
figure; 
boxplot(vertcat(dat.val),vertcat(dat.reg));
xlabel('Concentration / Region');
ylabel('Raw Intensity');

figure; scatter(opts.conc,vertcat(dat.sum),80,'red','o','filled'); 
xlim([min(opts.conc)-eps max(opts.conc)+eps]);
xlabel('Concentration / Region');
ylabel('Summed Intensity');

% Transform the roi into a decent alpha map...
alp = sum(roi,3) / max(roi(:));
idx = alp == 0;
alp(idx) = 0.25;

% Update the image for the alpha map
set(fig,'AlphaData',alp);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin)
% Read the arguments and then the data if it wasn't passed

% These are the defaults
opts.fP     = '/Users/jmckenzi/Desktop/';
opts.fN     = 'testRow21.imzML';
opts.mz     = [200 1000];
opts.numR   = 5;
opts.conc   = [];


% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    if strcmpi('mz',argsin{i})
        opts.mz = argsin{i+1};
    
    elseif strcmpi('conc',argsin{i})
        opts.conc = argsin{i+1};
        
    elseif strcmpi('regions',argsin{i})
        opts.numR = argsin{i+1};
        
    elseif strcmpi('folder',argsin{i})
        opts.fP = argsin{i+1};

    elseif strcmpi('file',argsin{i})
        opts.fN = argsin{i+1};
    
    end
end

if isempty(opts.conc)
    opts.conc   = 1:opts.numR;
end

if numel(opts.conc) ~= opts.numR
    error('Concentrations and number of regions are not the same size');
end

if ~exist([opts.fP opts.fN],'file')
    error('File does not exist');
end

if numel(opts.mz) ~= 2
    error('Give two mz values over a range');
end

if opts.mz(2) < opts.mz(1)
    error('Specify mz range properly');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz2] = getMZvectors(fP,fN,opts)
% Open the imzML files and get the mz vectors of each of the pixels

% imzML converter required
javaclasspath(deSlash([pwd ...
    '/packages/imzMLConverter/imzMLConverter.jar']));

% Structure to store all the mz values
mz = struct('qty',[],'vals',[],'inds',[]);
    
% File name
tmp = [fP fN];
disp(tmp);

imzML = imzMLConverter.ImzMLHandler.parseimzML(tmp);

nColumns = imzML.getWidth();
nRows    = imzML.getHeight();

% Temp place to store this file's mz values
%mzTmp = cell (nColumns*nRows-1,1);
%mzInd = zeros(nColumns*nRows-1,2);
%mzSze = zeros(nColumns*nRows-1,1);

% Make a grid just to plonk everything in
grid = zeros(nRows,nColumns);

% Run through the pixels...
i = 0;
for y = nRows-1:-1:1
    for x = 1:nColumns

        % Increase counter
        i = i + 1;     
        
        % Try, or skip
        try
            tmp = imzML.getSpectrum(x,y);
        catch err
            err
            continue;
        end

        % Skip empty scans
        if isempty(tmp)
            continue; 
        end 

        % Get the data
        xtmp       = imzML.getSpectrum(x,y).getmzArray();
        ytmp       = imzML.getSpectrum(x,y).getIntensityArray();

        [fx] = find(xtmp > opts.mz(1) & xtmp < opts.mz(2));

        if numel(fx) > 0
            %mzTmp{i,1} = [xtmp(fx) ytmp(fx)];
            %mzInd(i,:) = [x y];
            %mzSze(i,1) = size(mzTmp{i,1},1);     
            
            grid(y,x) = sum(ytmp(fx));
        end

    end
end

% Set first row to be the minimum value...
minVal = min(grid(:));
grid(1,:) = minVal;

mz.qty  = i;
%mz.vals = mzTmp;
%mz.inds = mzInd;
%mz.size = mzSze;
   
% Now that all of the files have been read, need to concatenate the various
% vectors into a format that is suitable for the mspmatch function, i.e. a
% n x 1 cell array
%mz2.allMz = vertcat(mz.vals);
%mz2.allIn = vertcat(mz.inds);
%mz2.allSz = vertcat(mz.size);
mz2.grid  = grid;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz] = getFromRawFile(fP,fN,opts)

% Structure to store all the mz values
%mz = struct('qty',[],'vals',[],'inds',[]);
    
% File name
tmp = [fP fN];
disp(tmp);

% Read in the entire data...
[sp,~,~,xy2D] = desiReadRaw([fP fN],false);

% Convert it to an image
[data] = desiReadRaw2Image(sp,xy2D);

% Now we need to extrac the m/z values we want to see
img = zeros(size(data));
for i = 1:size(data,1)
    for j = 1:size(data,2)
        
        % Temp data
        tmp = data{i,j};
        
        % mz mask
        mask = tmp(:,1) > opts.mz(1) & tmp(:,1) < opts.mz(2);
        
        % Sum region into img(i,j)
        img(i,j) = sum(tmp(mask,2));
        
    end
end

mz.qty = numel(data);
mz.grid = img;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dat2txt(fP,fN,dat,opts)
% Write out the data to a text file...

% Generate the file name...
px = strfind(fN,'.');
pN = [fN(1:px-1) '-Calib.txt'];

numF = size(dat,2);

maxS = 0;
for n = 1:numF
    if numel(dat(n).val) > maxS
        maxS = numel(dat(n).val);
    end
end

all = NaN(maxS,numF);
for n = 1:numF
    all(1:numel(dat(n).val),n) = dat(n).val;
end

all = [opts.conc; all];
dlmwrite([fP pN],all,'delimiter','\t','precision','%.6f');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

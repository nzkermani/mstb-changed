function dpnIMZML(~,~,fig,defP)
%dpnIMZML - load an imzML file!

f0 = findobj('Tag','manipulation');
close(f0);
f0 = findobj('Name','annotation');
close(f0);

% Ask the user for a file
[file.nam,file.dir,flag] = uigetfile({'*.imzML; *.mat'},'File Select',defP);
if flag == 0
    return
else
    file.ext = fileExtension(file.nam);
end

% Clear the axes... but what if we have annotation boxes?
sz = size(get(fig.ax.opt(2),'CData'));
set(fig.ax.opt(2),'CData',ones(sz));
sz = size(get(fig.ax.ms1(2),'CData'));
set(fig.ax.ms1(2),'CData',ones(sz));
sz = size(get(fig.ax.ms2(2),'CData'));
set(fig.ax.ms2(2),'CData',ones(sz));
sz = size(get(fig.ax.fu(2),'CData'));
set(fig.ax.fu(2),'CData',ones(sz));

f0 = findobj('Type','patch');
delete(f0);

% Now so let's read in the file
disp([file.dir file.nam]);

if strcmp(file.ext,'mat')
    % Then this is the matlab loading section.  We need to write a function
    % for it to be loaded back to the window.
    dpnMatlabLoad(file,fig);
    return
end

% Get the handle and then size of the imzML file
imzML = imzMLConverter.ImzMLHandler.parseimzML([file.dir file.nam]);

nCol = imzML.getWidth();
nRow = imzML.getHeight();

% These are the options that we will use at the beginning
[opts,flag] = getOptions(file);
opts.interp = 'linear';
if ~flag
    return;
end

numPix = nCol * (nRow-1);

sp.mz       = cell(numPix,1);%nCol,nRow-1);
sp.counts   = cell(numPix,1);%nCol,nRow-1);
sp.mzcounts = zeros(nCol,nRow-1);
sp.X        = zeros(numPix,1);%nCol,nRow-1);
sp.pixel    = zeros(numPix,2);

i = 0;

% Another waitbar
wb = waitbar(0,'Importing imzML');

% Loop through each pixel
for y = 1:nRow-1%nRow-1:-1:1
    for x = 1:nCol
        
        % Increase counter
        i = i + 1;
            
        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end 
            
        % Get the data
        xtmp = imzML.getSpectrum(x,y).getmzArray();
        ytmp = imzML.getSpectrum(x,y).getIntensityArray();
                
        % Round off m/z values
        sp.mz{i,1} = round(xtmp ./ opts.mzRes);        

        % Remove duplicated m/z values after rounding off
        [sp.mz{i,1},sp.counts{i,1}] = remDuplmz(sp.mz{i,1},ytmp);
        
        % Save the information...
        sp.mzcounts(x,y) = length(sp.counts{i,1});

        % Pixel location
        sp.pixel(i,:) = [x y];        
        
        % Count information
        %sp.counts{x,y} = ytmp(:,2);
        sp.X(i,1) = sum(log2(sp.counts{i,1}(sp.mz{i,1} > 500 ) + 20));
         
    end
    waitbar(y/(nRow+5),wb,'Importing imzML');
end

waitbar((y+1)/(nRow+5),wb,'Processing');

% With the 'raw' intensities of unaligned data, determine which is which
% mode
[ion,splitFig,flag] = dpnDetIonMode(sp);
graphFormat([file.dir 'PosNegSplit'],'png');
if flag
    disp('File processing failed');
    assignin('base','sp',sp);
    delete(wb);
    return
end
close(splitFig);
sp.mode = repmat(ion',[nCol 1]);

% Now we split up the data into two matrices
sp.mzcounts = reshape(sp.mzcounts,[nCol*(nRow-1) 1]);
[sp1,sp2] = splitPosNeg(sp);

% Now we run the m/z alignment function on the two sp matrices separately.
% From now on they are entirely separate datasets which will have different
% m/z vectors (intensity is a different matter!)
[MZ1,X1] = finalProc(sp1,nRow,nCol,opts);
[MZ2,X2] = finalProc(sp2,nRow,nCol,opts);

% Once you've imported the data, we should look at interpolating in the
% empty rows
[X1,intM1,ii1] = dpnPerformInterpolation(X1,opts.interp,[]);
[X2,intM2,ii2] = dpnPerformInterpolation(X2,intM1,[]);

% Create a guidata structure, and then save the results, and then display
% them in the right boxes...
dpn.file = file;
dpn.opts = opts;
dpn.opts.interp = intM2;

dpn.d1.mz = MZ1;
dpn.d1.sp = X1;
dpn.d1.isInterp = ii1;

dpn.d2.mz = MZ2;
dpn.d2.sp = X2;
dpn.d2.isInterp = ii2;

dpn.fig = fig;
dpn.defP = defP;

dpn.mode = 'dual';

delete(wb);

% Add the guidata...
guidata(fig.fig,dpn);

% Update the images...
dpnUpdateMS([],[],fig,'force');

% Update the fused ion image to end...
dpnUpdateFusedIonImage([],[],fig);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts,flag] = getOptions(file)

% Decide whether we can continue; if false then we should abort
flag = true;

% Options for processing imzML files 
opts.mzFrac   = 0.01;
opts.mzRes    = 0.002;
opts.ppmRes   = 8;
opts.mzRange  = [100 1000];
opts.prompt = {'Fraction of empty pixels per m/z',...
    'm/z resolution',...
    'Maximum ±ppm shift',...
    'm/z range'};

% These are the defaults
defs = {num2str(opts.mzFrac),...
    num2str(opts.mzRes),...
    int2str(opts.ppmRes),...
    num2str(opts.mzRange)};
        
% Ask the user the question
resp = inputdlg(opts.prompt,'Processing parameters',1,defs);

if isempty(resp)
    % User aborted processing
    disp('User aborted the processing');
    opts = [];
    flag = false;
    return
end
    
% Convert to numbers - no error checking/validation though
opts.mzFrac   = str2double(resp{1});
opts.mzRes    = str2double(resp{2});
opts.ppmRes   = str2double(resp{3});            
opts.mzRange  = str2num(resp{4});

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,val] = remDuplmz(mz,val)

diffmz           = diff(mz)';
[~, mzindcs] = find(diffmz>1);
dubmzindcs       = find(diff(mzindcs)>1);

if ~isempty(dubmzindcs)
    counts_temp = val([1 mzindcs+1]);
    for i = dubmzindcs
        counts_temp(i+1) = sum(val(mzindcs(i)+1:mzindcs(i+1)));
    end
    mz  = mz([1 mzindcs+1]);
    val = counts_temp;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp1,sp2] = splitPosNeg(sp)
% Split the things into the two modes

sz = numel(sp.mode);
mode = reshape(sp.mode,[sz 1]);

% Indices of pixels defining the split
m1 = mode == 1;
m2 = mode == 2;

% Fieldnames
%fn = fieldnames(sp);
fn = {'mz','counts','mzcounts','X','pixel'};

for n = 1:numel(fn)
    
    % Get the bit of it
    tmp = sp.(fn{n});
    
    sp1.(fn{n}) = tmp(m1,:);
    sp2.(fn{n}) = tmp(m2,:);
    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MZ,X] = finalProc(sp,nRow,nCol,opts)

% Reshape the matrices so now we just have a long list... then extract them
% into the two modes
%sp.objectmz     = reshape(sp.mz,1,(nRow-1)*nCol);
%sp.X            = reshape(sp.X,(nRow-1)*nCol,1);
%sp.mode         = reshape(sp.mode,(nRow-1)*nCol,1);
%sp.mzcounts     = reshape(sp.mzcounts,(nRow-1)*nCol,1);

% mzcounts is just the number of data points per spectrum
sp.objmzcounts = sp.mzcounts;

% Cumulative sum of vector lengths
mzobjindcs = cumsum([1 sp.objmzcounts']);
 
% Determine the mz values throughout the sample
objectmzs = zeros(1,sum(sp.objmzcounts));
nObjPixels = sum(sp.X > 1);
for i = 1:nObjPixels
    objectmzs(mzobjindcs(i):mzobjindcs(i+1)-1) = sp.mz{i,1}; 
end

% Find all unique values of the mz vector, then the frequency of each
% throughout the sample
MZ = unique(objectmzs); 
n  = hist(objectmzs,MZ); 
clear objectmzs;

% Filter out mz values that don't appear very frequently
MZ = MZ(n > nCol*(nRow-1) * opts.mzFrac); 
X  = zeros(nRow-1,nCol,length(MZ));

for i = 1:size(sp.pixel,1)
        
        % What is the (x,y) location of this pixel?
        x = sp.pixel(i,1);
        y = sp.pixel(i,2);
        
        % This pixel's mz vector
        mz  = sp.mz{i,1}; 
        
        % Align mz intensities to a common MZ feature vector
        [~,mzindcs,MZindcs] = intersect(mz,MZ); 
        
        % Put into the new matrix
        X(y,x,MZindcs)    = sp.counts{i,1}(mzindcs);

end

% Restore the mz scale to a meaningful set of numbers
MZ = MZ * opts.mzRes;

% Combine peaks that were split
[X,MZ] = combSplitPeaks2(X,MZ,opts.ppmRes);

% Scale median value to 40 for consistency
%medX = median(X(X~=0));
%X = 40 * X ./ medX; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



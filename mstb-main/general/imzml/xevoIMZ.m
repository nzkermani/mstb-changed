function [ op ] = xevoIMZ(varargin)
% imzMLreader - read in the data from an imzML file.
%
% How to use this marvellous function?
% Type the following to make it work, you'll need to change the things that
% appear between <> 
%
% imzMLreader('folder',<a>,'file',<b>,'mz',<c>);
%
% <a>, folder where the file is
% <b>, name of the file (imzML)
% <c>, m/z range to be imported, e.g. [260.16 260.17]


% Get the addition stuff
[opts] = readArgsData(varargin);

% Get the stuff from the file
[op] = getMZvectors(opts.fP,opts.fN,opts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = readArgsData(argsin)
% Read the arguments and then the data if it wasn't passed

% These are the defaults
opts.fP     = '/Users/jmckenzi/DB/Herbs/imzML/';
opts.fN     = '2015-09-08_Thyme_Neg.imzML';
opts.inst   = 'xevo';

%opts.fP     = '/Users/jmckenzi/DB/3D-Liver/ProfileData/';
%opts.fN     = '7TopRight, 17BottomRight, 27BottomLeft, 37TopLeft-profile.imzML';
%opts.inst   = 'orbi';

opts.mz     = [600 1000];

% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    if strcmpi('mz',argsin{i})
        opts.mz = argsin{i+1};
        
    elseif strcmpi('folder',argsin{i})
        opts.fP = argsin{i+1};

    elseif strcmpi('file',argsin{i})
        opts.fN = argsin{i+1};
    
    end
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
function [new] = getMZvectors(fP,fN,opts)
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
mzTmp = cell (nColumns*nRows-1,1);
mzInd = zeros(nColumns*nRows-1,2);
mzSze = zeros(nColumns*nRows-1,1);
%mzPPM = cell (nColumns*nRows-1,1);

% Make a grid just to plonk everything in
grid = zeros(nRows,nColumns);

% Run through the pixels...
i = 0;
for y = nRows-1:-1:1
    for x = 1:nColumns


        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end 

        % Get the data
        xtmp       = imzML.getSpectrum(x,y).getmzArray();
        ytmp       = imzML.getSpectrum(x,y).getIntensityArray();

        [fx] = find(xtmp > opts.mz(1) & xtmp < opts.mz(2));% & ytmp > 0);
        
        % Store the vectors in a cell array
        if numel(fx) > 0        
            
            % Increase counter
            i = i + 1;            

            mzTmp{i,1} = [xtmp(fx) ytmp(fx)];
            mzInd(i,:) = [x y];
            mzSze(i,1) = size(mzTmp{i,1},1);     

            % Total ion image
            grid(y,x) = sum(ytmp(fx));
        else
            fx
        end

    end
end

mz.qty  = i;
mz.vals = mzTmp(1:i);
mz.inds = mzInd(1:i,:);
mz.size = mzSze(1:i,:);
clear mzTmp mzInd mzSze;

% Let's process the data in full profile mode, so no variable extraction
% here.  Thankfully, it appears that all m/z values across the sample are
% somewhat quantized (hopefully true for all files) so we can just find all
% unique m/z values and create the full data matrix accordingly

% Concat all vectors together to get the m/z vector
allMZ = vertcat(mz.vals{:});

% Orbitrap and Xevo mz values are different; Orbi produces unique values
% but Xevo gives a series of unique values in spite of being recorded at
% seemingly infinite precision
switch opts.inst
    case 'xevo'
        unqMZ = unique(allMZ(:,1));
        numMZ = numel(unqMZ);
        clear allMZ;
        
    case 'orbi'
        % Here we calculate the ion differences between datapoints
        [fit] = detMZdiffs(mz.vals,false);
        
        unqMZ = varMZvector(opts.mz(1),opts.mz(2),fit);
        numMZ = numel(unqMZ);
end

% Create a data matrix of that size...
img = zeros(nRows,nColumns,numMZ);
size(img)

% Now place the data into the matrix
tic
switch opts.inst    
    % XEVO: this uses the intersect method
    case 'xevo'        
        for n = 1:mz.qty
            % Use intersect to find matching values...
            [~,orig,new] = intersect(mz.vals{n}(:,1),unqMZ);
            
            % Add the intensities into the mega matrix
            img(mz.inds(n,2),mz.inds(n,1),new) = mz.vals{n}(orig,2);
        end
        
    case 'orbi'
        for n = 1:mz.qty
            % Use the interpolation method to perform this task. First 
            % convert 0 -> NaN as it buggers up the interpolation at the
            % edges if there is a big drop/rise in intensity
            vals = mz.vals{n}(:,2);
            %vals(vals == 0) = NaN;            
            
            % Begin the interpolation
            ival = interp1(mz.vals{n}(:,1),mz.vals{n}(:,2),unqMZ);
            %ival(isnan(ival)) = 0;
            
            % Add in the right place
            img(mz.inds(n,2),mz.inds(n,1),:) = ival;
            
            if mod(n,50) == 0
                disp([int2str(n) '/' int2str(mz.qty)]);
            end
        end
end
toc

assignin('base','img',img);
assignin('base','unqMZ',unqMZ);
    
% From here we will do the neighbouring correlation function that will 
% remove peaks from the matrix
[new.mz,new.img,new.opts] = corrIMZ(unqMZ,img);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

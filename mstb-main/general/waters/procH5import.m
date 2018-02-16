function [ gc ] = procH5import(file)
%procH5import - load all of the files into memory / or sequentially for
%stuff if required

% File name
if nargin == 0
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/H5-Test-170510-135411.h5';
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/Windows-H5-Test-170511-101632-SubFile.h5';
    
    file = 'E:\Box Sync\H5 Proc\DESI-Raman-170525-113747.h5';
    %file = 'E:\Box Sync\H5 Proc\CRUK-GC-Neg-Final170602-100814.h5';
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/PP-FullRange-10-Proc-170522-154228.h5';
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/Windows-H5-Test-170511-101632.h5';
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/Mac-IMZML-Test-170524-143758.h5';
    %file = '/Users/jmckenzi/Documents/Box Sync/H5 Proc/Mac-H5-Test-PEAK-THRESH-2-170527-020429.h5';
end

% Work out how many datasets are included within...
info = h5info(file);
numF = size(info.Groups,1);

% Read in the mz vector
gc.fn = cell(numF,1);
gc.mz = h5read(file,'/mz');

% Read in each part...
sp = cell(numF,4);
for n = 1:numF
    
    tic
    gc.fn{n} = h5readatt(file,['/file' int2str(n)],'fileName');
    tmp = h5read(file,['/file' int2str(n) '/img']);
    rawMZ = h5read(file,['/file' int2str(n) '/rawMZ']);
    
    %mzImg = h5read(file,['/file' int2str(n) '/mzImg']);
    %assignin('base','mzImg',mzImg);
    
    sz = size(tmp);
    
    % Reshape it a little
    tmp = reshape(tmp,[sz(1)*sz(2) sz(3)]);
    
    % Determine the background here...
    %[tobg] = procH5background(gc.mz,tmp,sz);
    tobg = true(size(tmp,1),1);
        
    % Create an index to determine which file the pixels belong to
    id = ones(sz(1)*sz(2),1) * n;
    
    % Save in sp
    sp{n,1} = tmp;
    sp{n,2} = sz;
    sp{n,3} = id;
    sp{n,4} = tobg;
    
    toc
    disp([int2str(n) '/' int2str(numF)]);
    
end

% Format into a single matrix rather than the cell array.
gc.sp = vertcat(sp{:,1});
gc.sz = vertcat(sp{:,2});
gc.idx = vertcat(sp{:,3});
gc.tobg = vertcat(sp{:,4});

end


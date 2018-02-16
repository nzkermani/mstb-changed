function [vals] = imzmlp(hits)
% Do some parsing of imzml files in a specified folder & subfolder. Use the
% function fileFileType to get the imzML files and locations

% imzML converter required
javaclasspath(['/Users/jmckenzi/Github/jsmCode/'...
    'imzMLConverter/imzMLConverter.jar']);

% How many files are there?
numF = size(hits,1);

vals = cell(numF,4);

% Loop through the files...
for n = numF-10:numF
    
    % File name
    tmp = [hits{n,1} filesep hits{n,2}];
    %disp([fN{n} '-centroid.imzML']);
    
    %tic
    
    [vals{n,2},vals{n,3},vals{n,4}] = parseXML(tmp);
    vals{n,1} = hits{n,2};
    
end


end


function [] = abc()
    imzML = imzMLConverter.ImzMLHandler.parseimzML(tmp);
    
    nColumns = imzML.getWidth();
    nRows    = imzML.getHeight();
    
    % Temp place to store this file's mz values
    mzTmp = cell (nColumns*nRows-1,1);
    mzInd = zeros(nColumns*nRows-1,3);
    mzSze = zeros(nColumns*nRows-1,1);

    % Run through the pixels...
    i = 0;
    for y = nRows-1:-1:1
        for x = 1:nColumns
            
            % Increase counter
            i = i + 1;            
            
            % Skip empty scans
            if isempty(imzML.getSpectrum(x,y))
                continue; 
            end 
            
            % Get the data
            xtmp       = imzML.getSpectrum(x,y).getmzArray();
            ytmp       = imzML.getSpectrum(x,y).getIntensityArray();
            
            [fx] = find(xtmp > mzVal(1) & xtmp < mzVal(2));
            
            if numel(fx) > 0
                mzTmp{i,1} = [xtmp(fx) ytmp(fx)];
                mzInd(i,:) = [n x y];
                mzSze(i,1) = size(mzTmp{i,1},1);                
            end

        end
    end
    %mz(n).qty  = i;
    %mz(n).vals = mzTmp;
    %mz(n).inds = mzInd;
    %mz(n).size = mzSze;

    toc
    
    % Want to export the raw data into the h5 file. However, h5 files don't
    % support cell array/structures, so we need to convert all the data
    % into a single double-like matrix
    varName = ['desi/raw/s' int2str(n) '/'];
    
    % Concatenate the mz values and intensities
    tmpSpec = vertcat(mzTmp{:});
    h5W(h5P,[varName 'spec'],tmpSpec,1000,2);
    
    % Write the sizes of each pixel's spectrum
    h5W(h5P,[varName 'size'],mzSze,1000,2);
    
    % Write the pixel information of the spectra
    h5W(h5P,[varName 'info'],mzInd,1000,2);
    
    % Attributes: file name and path
    h5writeatt(h5P,['/' varName],'imzML-path',fP);
    h5writeatt(h5P,['/' varName],'imzML-name',[fN{n} '-centroid.imzML']);
        
    % Here we format the data to give it its local mz vector (lmz), and
    % then we do the diagnostics to have a look...
    %     [MZ,X] = micheleImport(mzTmp,mzInd,...
    %         'Method','tobg',...
    %         'combSplit','fixed',...
    %         'mzRes',0.001);
    %     diagImport(mzTmp,mzInd,MZ,X,[false false true]);
    
    % Here we format the data using a variable ppm based combination method
    % for the split peaks
    [MZ,X] = micheleImport(mzTmp,mzInd,...
        'Method','tobg',...
        'mzRes',0.001,...
        'combSplit','ppm',...        
        'ppmRes',5);
    
    %diagImport(mzTmp,mzInd,MZ,X,[false true false]);
    %default params were false true true
    
    % Now write these two things to the h5 file as well. Also need to add
    % any of the processing parameters that have been used...
    varName = ['desi/lmz/s' int2str(n) '/'];    
    h5W(h5P,[varName 'x'], X, 1000,2);
    h5W(h5P,[varName 'mz'],MZ,1000,2);
    
    % Finally, what about the optical image? Whilst we are here we should
    % certainly read them in...
    tmp = [fP fN{n} filesep fI{n}];
    disp(fI{n});
    
    % Import image
    opImg = imread(tmp);
    
    % Save to the h5 file
    h5W(h5P,['opt/s' int2str(n)],opImg,1000,2);
    
    % Attributes to be written
    h5writeatt(h5P,['/opt/s' int2str(n)],'opt-path',[fP fN{n} filesep]);
    h5writeatt(h5P,['/opt/s' int2str(n)],'opt-name',fI{n});

    


% So we shouldn't need to return anything from this function except the
% processing parameters used in the above stages

%return

% Now that all of the files have been read, need to concatenate the various
% vectors into a format that is suitable for the mspmatch function, i.e. a
% n x 1 cell array
mz2.allMz = vertcat(mz.vals);
mz2.allIn = vertcat(mz.inds);
mz2.allSz = vertcat(mz.size);
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


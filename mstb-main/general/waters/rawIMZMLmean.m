function [data,meanSpec,meanFreq] = rawIMZMLmean(file,mzVec)
%rawH5mean - determine mean spectrum

% What about the ppm method...
if mzVec(2)-mzVec(1) == mzVec(3)-mzVec(2)
    method = 'fixed';
    mzFunction = @methodFixed;
    res = mzVec(2)-mzVec(1);
else
    method = 'variable';
    mzFunction = @methodVariable;
    res = [];
end

wb = waitbar(0,'Hello');

% Add the package path
javaclasspath('packages/imzMLConverter/imzMLConverter.jar');

% Get the handle and then size of the imzML file
imzML = imzMLConverter.ImzMLHandler.parseimzML(file);
nCol = imzML.getWidth();
nRow = imzML.getHeight();

i = 0;

% Save all the data : x|y|mz,sp
data = cell((nRow-1) * nCol,3);
meanSpec = zeros(size(mzVec));
meanFreq = zeros(size(mzVec));

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
        mz = imzML.getSpectrum(x,y).getmzArray();
        sp = imzML.getSpectrum(x,y).getIntensityArray();
        
        % mz mask
        mask = mz > mzVec(1) & mz < mzVec(end);

        % Perform matching
        [a] = mzFunction(mzVec,mz(mask),res);
            
        % Add in to the matrix
        meanSpec(a,1) = meanSpec(a,1) + sp(mask);
        meanFreq(a,1) = meanFreq(a,1) + 1;
        
        % Do this as well
        data{i,1} = x;
        data{i,2} = y;
        data{i,3} = [mz sp];
        
    end
    waitbar(y/nRow,wb);
end

delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a] = methodFixed(mzVec,mzSp,res)
% Fixed mode binning, usually at high resolution

a = round((mzSp-mzVec(1)) / res) + 1;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a] = methodVariable(mzVec,mzSp,res)
        
[a,~] = comp2vec(mzSp',mzVec');
            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,Sp] = getMSImageCentroidMatrix(imzML,mzRange)
% getMSImageDataMatrix exctracts MS image intensity matrix from..
% imzML data object storing a data acquired in centroid mode 
%    Input:  imzML  - the MS image data object of a specimen 

%    Output: MZ - mz values 
%            X  - output matrix [rows x columns x mz]
% Author: Kirill A. Veselkov, Imperial College London 2012. Modified by
% Nazanin Z. Kermani, , Imperial College London 2016.

% mzRange is specified as second input argument.
if nargin == 1
    mzRange = [0 Inf];
end
if isempty(mzRange)
    mzRange = [0 Inf];
end

nColumns = imzML.getWidth();
nRows    = imzML.getHeight();

Sp = zeros(nRows-1,nColumns,30000);
mz = zeros(nRows-1,nColumns,30000);

h = waitbar(0,'Reading imzML file....');
steps = (nRows-1)*nColumns;
for y = nRows-1:-1:1
    for x = 1:nColumns
        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end
        
        % Get the data
        imz        = imzML.getSpectrum(x,y).getmzArray();
        counts     = imzML.getSpectrum(x,y).getIntensityArray();
        
        fx = imz >= mzRange(1) & imz <= mzRange(2);
        
        lengthPks  = sum(fx);%length(pks);
        Sp(y,x,1:lengthPks) = counts(fx);
        mz(y,x,1:lengthPks) = imz(fx);
            waitbar((steps-(((y-1)*(nRows-1))+y)) / steps);
    end
end
close(h);
Sp = Sp(:,:,any(mz,3));
mz = mz(:,:,any(mz,3));

end


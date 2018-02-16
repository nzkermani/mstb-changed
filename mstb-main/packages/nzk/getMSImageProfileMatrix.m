function [mz,Sp] = getMSImageProfileMatrix(imzML,options)
%% getMSImageDataMatrix exctracts MS image intensity matrix from..
%% imzML data object storing a data acquired in centroid mode 
%    Input:  imzML  - the MS image data object of a specimen 

%    Output: MZ - mz values 
%            X  - output matrix [rows x columns x mz]
%% Author: Kirill A. Veselkov, Imperial College London 2012. Modified by
% Nazanin Z. Kermani, , Imperial College London 2016.

if nargin < 2
     options = [];
end

options = getVarArginMS(options); 

nColumns = imzML.getWidth();
nRows    = imzML.getHeight();

Sp = zeros(nRows-1,nColumns,30000);
mz = zeros(nRows-1,nColumns,30000);
for y = nRows-1:-1:1
    for x = 1:nColumns
        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end 
        % Get the data
        imz        = imzML.getSpectrum(x,y).getmzArray();
        counts     = imzML.getSpectrum(x,y).getIntensityArray();
        lengthPks  = numel(counts);%length(pks);
        Sp(y,x,1:lengthPks) = counts;
        mz(y,x,1:lengthPks) = imz;
    end
end

Sp = Sp(:,:,any(mz,3));
mz = mz(:,:,any(mz,3));

end

%% 
function options = getVarArginMS(argsin)
options.display   = 0;
options.mzrange   = [200 1000];

nArgs = length(argsin);
for i=1:2:nArgs
    if strcmp('display',argsin{i})
        options.display   = argsin{i+1};
    elseif strcmp('mzrange',argsin{i})
        options.mzrange   = argsin{i+1};
    end
end

end
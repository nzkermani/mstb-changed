function [mz,Sp] = getMSImageProfileMatrix(imzML,options)
%% getMSImageDataMatrix exctracts MS image intensity matrix from..
%% imzML data object with a given resolution (default 0.001).
%    Input:  imzML  - the MS image data object of a specimen 

%    Output: MZ - mz values 
%            X  - output matrix [rows x columns x mz]
%% Author: Kirill A. Veselkov, Imperial College London 2012. 
warning off all;   currentFolder = pwd;
cd(currentFolder); addpath(genpath(currentFolder));

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
        imz         = imzML.getSpectrum(x, y).getmzArray();
        counts     = imzML.getSpectrum(x, y).getIntensityArray();
        [pks,locs] = findpeaks(counts);
        imz        = imz(locs);
        Sp = [imz'; pks'];
        fname = sprintf('ms%d%d.txt', y, x);
        fileID = fopen(fname,'w');
        fprintf(fileID,'%12.8f \t %12.8f\n',Sp);
        fclose(fileID);
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
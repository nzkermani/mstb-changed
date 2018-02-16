function [ output_args ] = readWatersHDF(file)
%readWatersHDF - import the data from an HDF file created from a Waters Raw
%file. Hopefully it will be quicker than using the Waters library each
%time...

% Read in the easy bits, such as TIC / totPts / xy2D
totInt = h5read(file,'/totInt');
totPts = h5read(file,'/totPts');
xy2D = h5read(file,'/xy2D');

% Determine an average spectrum somehow?
mzRange = [50 1500];
res = 0.01;
[sp] = avgSpec(file,numel(totInt),mzRange,res);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,sp] = avgSpec(file,numS,mzRange,res)

% Add a simple buffer
mzRange = mzRange + [-10 10];

mz = mzRange(1):res:mzRange(2);
numM = numel(mz);
sp = zeros(numM,1);

for n = 1:numS
    
    % Read in the scan
    tmp = h5read(file,['/raw/scan/' int2str(n)]);
    
    % Convert mz values to indices
    idx = floor(tmp(:,1) / res) - (mz(1) / res) + 1;
    
    % Add to the spectrum
    sp = sp + accumarray(idx,tmp(:,2),[numM 1]);
   
end

% Trim out zeros from the beginning. Want to have first intensity being non
% zero...
fx = find(sp > 0,1,'first');
fy = find(sp > 0,1,'last');

mz = mz(fx:fy);
sp = sp(fx:fy);




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mznew,spnew ] = watersAnalyteFile(file,data)
%watersAnalyteFile - import the txt analyte file

% Create a waitbar
wb = waitbar(0.25,'Importing Analyte File');

% Read in unless provided
if ~isempty(file) && isempty(data)    
    data = dlmread(file);    
end

% Change
waitbar(0.75,wb,'Processing Analyte File');

% Extract the various bits from the file
x  = data(4:end,2);
y  = data(4:end,3);
mz = data(3,:)';
sp = data(4:end,4:end);

% Trim dodgy variables
mask = mz > 0;
mz = mz(mask);
sp = sp(:,mask);

% Need to convert x/y into a regular grid...
xy2D = get2Dcoord(x,y);
[nRows,nCols] = size(xy2D);
xy2D = reshape(xy2D,[nRows * nCols 1]);

% How about sorting m/z values
[mznew,idx] = sort(mz);

% Reorder...
%spnew = zeros(size(sp,1),numel(mz));
spnew = sp(xy2D,idx);
spnew = reshape(spnew,[nRows nCols numel(mz)]);
spnew = flipud(spnew);

% Delete waitbar
delete(wb);

end


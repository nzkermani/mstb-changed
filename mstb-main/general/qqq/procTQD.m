function [img] = procTQD(tqd)
%procTQD - process the TQD structure to make images!

% What are the dimensions of the images to be created?
dim1 = max(tqd.data(:,1))
dim2 = max(tqd.data(:,3))+1

% How many channels?
numC = numel(tqd.head)

% Define the x/y coordinates for the correct image orientation
if dim1 > dim2
    xc = tqd.data(:,3) + 1;
    yc = tqd.data(:,1);
    
    % Make an empty image
    img = NaN(dim2,dim1,numC);

else
    xc = tqd.data(:,1);
    yc = tqd.data(:,3) + 1;
    img = NaN(dim1,dim2,numC);
end

% Now run through the data and put all in the right places in the matrices
for n = 1:size(tqd.data,1)
    
    x = xc(n);
    y = yc(n);
    z = tqd.data(n,4:end);
    
    img(x,y,:) = z;
    
end

end


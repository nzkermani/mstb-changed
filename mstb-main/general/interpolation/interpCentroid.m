function [ output_args ] = interpCentroid( input_args )
%interpCentroid - interpolate centroid mode data to demonstrate how it
%fails with peak positional variation.

mz1 = 500.4:0.001:500.6;
mz2 = 500.4:0.01:500.6;

% These are the locations of the existing peaks...
peak = 500.495:0.0005:500.505;
numP = numel(peak);

for n = 1:numP
    
    
    newMZ = [peak(n)-0.01:0.005:peak(n)+0.01];
    newInt = zeros(size(newMZ));
    sz = (numel(newMZ)+1)/2
    newInt(sz) = 1;
    
    aa = interp1(newMZ, [0 1 0],mz1);
    
end

end


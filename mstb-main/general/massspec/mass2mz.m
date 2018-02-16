function [mzgrid] = mass2mz(mass,adl)
% Convert a series of mass values from a database into m/z values dependent
% on the various ion forms

% How many adducts are there?
numAdd = size(adl.list,1);
mzgrid = zeros(numel(mass),numAdd);

for r = 1:numAdd
    mzgrid(:,r) = (adl.list(r,3) * mass) + adl.list(r,1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


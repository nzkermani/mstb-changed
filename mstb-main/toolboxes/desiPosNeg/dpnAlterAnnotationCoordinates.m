function [dpn] = dpnAlterAnnotationCoordinates(dpn)
%dpnAlterAnnotationCoordinates - modify the coordinates a little

numA = size(dpn.anno,1)
for n = 1:numA
    
    % Check that it isn't a single pixel annotation
    sz = numel(dpn.anno{n,6})
    if sz == 1
        continue;
    end
    
    % Check if there are halves...
    dpn.anno{n,6} = [ceil(dpn.anno{n,6}(1)) floor(dpn.anno{n,6}(2))];
    dpn.anno{n,7} = [ceil(dpn.anno{n,7}(1)) floor(dpn.anno{n,7}(2))];

end


end


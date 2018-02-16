function [ output_args ] = dpnPaperResiduals(op)
%dpnPaperResiduals - for all ions, determine the residual between measured
%and interpolated intensities, and express as a ratio of the measured
%intensity.  In this way, we can generate a rough idea of where the
%residuals are largest, i.e. which pixels show the biggest deviation.

res = NaN(size(op.n1));

sz = size(op.crr)

for i = 1:sz(1)
    
    for j = 1:sz(2)
        
        % X and Y intensities
        x = squeeze(op.n1(i,j,:));
        y = squeeze(op.n2(i,j,:));
        
        % Set x to be measured and y interpolated
        if op.interp1(i,j)
            z = x;
            x = y;
            y = z;
        end
        
        % Determine the residual by subtracting the two sets...
        res(i,j,:) = x - y;
                
    end
    
end


res2 = reshape(res,[sz(1)*sz(2) size(res,3)]);


end


function [x]=ticNorm(x, method)
%% ticNorm Total Ion Current (TIC) normalization the intensity matrix 
%%    Input:  x  - Intensity matrix [rows x columns x mz]

%    Output: x  - Intensity matrix normalized [rows x columns x mz]
%      
%% Author: Nazanin Z. Kermani, , Imperial College London 2016.

if(method == 'TIC')
    if (ndims(x)>2)
        temp = reshape(x, size(x,1)*size(x,2), size(x,3));
        tic = sum(temp,2);
        temp = bsxfun(@rdivide,temp,tic);
        x = reshape(temp, size(x,1),size(x,2), size(x,3));
    end
    
end
end
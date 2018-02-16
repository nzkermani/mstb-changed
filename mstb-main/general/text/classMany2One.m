function [ newY ] = classMany2One(y,sep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numC = size(y,2);

numO = size(y,1);

newY = cell(numO,1);

if nargin == 1
    sep = '@';
end

for n = 1:numO
    
    tmp = y{n,1};
    
    if isnumeric(tmp)
        if isnan(tmp)
            tmp = 'NaN';
        else
            tmp = int2str(tmp);
        end
    end
    
    for r = 2:numC
        
        tmp2 = y{n,r};
        if isnumeric(tmp2)
            if isnan(tmp2)
                tmp2 = 'NaN';
            else
                tmp2 = int2str(tmp2);
            end
        end
        
        try
            tmp = [tmp{:} sep tmp2];
        catch
            tmp = [tmp sep tmp2];
        end
        
    end
    
    newY{n,1} = tmp;
    
    
end


end


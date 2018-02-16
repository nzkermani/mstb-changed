function [ fx ] = strfindMultiple(la,lb)
%strfindMultiple - find items in list A (la) that occur in list B (lb), note
%that this is a multi-way extension from strcmp

numF = numel(lb);

fx = false(numel(la),1);

for n = 1:numF
    
    tmp = strfind(la,lb{n});
    
    tmp = ~cellfun(@isempty,tmp);
    
    fx = fx + tmp;
    
end

fx = fx > 0;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [txt,all] = uniqueRows(anno)

numB = size(anno,2);

all = zeros(size(anno,1),numB);

for n = 1:numB
    
    [~,idx,tmp] = unique(anno(:,n));
    
    for r = 1:numel(idx)
        fx = tmp == r;
        all(fx,n) = idx(r);
    end    
end

% These are the unique numbers
unq = unique(all,'rows');

% Now convert these to text
txt = cell(size(unq));
for n = 1:size(unq,1)    
    for r = 1:size(unq,2)        
        txt(n,r) = anno(unq(n,r),r);        
    end    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

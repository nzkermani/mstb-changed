function [op] = mtspAlignAnnotations(data)
%mtspAlignAnnotations - align data according to annotations...


% Concatenate all annotations together...
allAnn = vertcat(data.annos);

% Combine together
comb = classMany2One(allAnn,'@');

% Determine uniques
[unq,~,ind] = unique(comb);
numA = numel(unq);
fileID = cell(size(data,2),1);
foldID = cell(size(data,2),1);

% Determine the frequency of annotations
frq = hist(ind,1:numA);

% Loop
for n = 1:size(data,2)
    
    % Combine
    comb = classMany2One(data(n).annos,'@');
    
    if isempty(comb)
        continue;
    end
    
    % Text intersect...
    [~,x,y] = intersect(unq,comb);
    
    % Now we have the indices to make a larger aligned data matrix.
    data(n).al = zeros(size(data(n).sp,1),numA);
    data(n).al(:,x) = data(n).sp(:,y);
    
    fileID{n,1} = repmat({data(n).file},[size(data(n).sp,1) 1]);
    foldID{n,1} = repmat({data(n).subType},[size(data(n).sp,1) 1]);
    
end
    
% Combine into a single matrix
al = vertcat(data.al);

% Combine histologies
histID = vertcat(data.histID);

% Combine files...
ie = cellfun(@isempty,fileID);
fileID = vertcat(fileID{~ie,:});
foldID = vertcat(foldID{~ie,:});

% Split the adduct from the annotation
locAt = cellfun(@strfind,unq,repmat({'@'},[numA 1]));
annos = cell(numA,2);
for n = 1:numA
    annos{n,1} = unq{n,1}(1:locAt(n)-1);
    annos{n,2} = unq{n,1}(locAt(n)+1:end);
end

% REturn the op structure
op.mz = (1:numA)';
op.sp = al;
op.annos = annos;
op.freq = frq;
op.meta.fileID = fileID;
op.meta.histID = histID;
op.meta.foldID = foldID;

end


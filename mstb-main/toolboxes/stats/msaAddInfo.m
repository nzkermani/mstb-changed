function [ msa ] = msaAddInfo(msa,info,heads)
%msaAddInfo - add new metadata information to the msa.meta structure. THis
%is originally for Olivia's data to insert dates into the endometrial
%samples to check for a batch effect.

% How many entries?
sz = size(msa.meta.fileID,1)

% Create empty bits in msa.meta structure...
numF = size(info,2) - 1;
extra = cell(sz,numF);

% Scroll through
for n = 1:size(info,1)
    
    % Find entries in the fileID list
    fx = ~cellfun(@isempty,strfind(msa.meta.fileID,info{n,1}));
    
    % Skip if none found
    if sum(fx) == 0
        continue;
    end
    
    extra(fx,:) = repmat(info(n,2:end),[sum(fx) 1]);
    
end

% Add into the structure using provided heads
if nargin == 2
    heads = [];
end
if isempty(heads)
    heads = cell(1,numF);
    for n = 1:numF
        heads{n} = ['extra' int2str(n)];
    end
end

% Add 'em in
for n = 1:numF
    
    if isnumeric(extra{1,n})
        msa.meta.(heads{n}) = cell2mat(extra(:,n));
    else
        msa.meta.(heads{n}) = extra(:,n);
    end
    
end
   
msa.meta

end


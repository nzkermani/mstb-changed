function [anno] = xxxAnnotationConsolidate(anno)
%xxxAnnotationConsolidate - ensure that the annotation labels are applied
%consistently throughout the annotation table in DESI/DESI pos/neg

% We will use the second column of numbers to be the arbiter as to group
% names across the file

% What are the labels in this file?
labels = cell2mat(anno(:,2));
[unq,~,ind] = unique(labels);
numG = numel(unq);

% Loop through each...
for n = 1:numG
    
    fx = ind == n;
    
    % Does this have multiple annotation name?
    names = anno(fx,5);
    [unqN,~,indN] = unique(names);
    numN = numel(unqN);
    
    % Question?
    if numN > 1
        
        % Find the most common name
        freq = hist(indN,1:numN);
        [~,comm] = sort(freq,'descend');
        new = unqN{comm(1)};
        
        % Check that it isn't 'ID x' format...
        if strcmp(new(1:3),'ID ')
            new = unqN{comm(2)};
        end
        
        % Determine indices of not most common names
        %idx = ~strcmp(names,new);
        
        % Rename the files accordingly...
        anno(fx,5) = {new};
        
    end
    
end

end


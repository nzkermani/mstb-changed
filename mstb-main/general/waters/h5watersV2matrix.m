function [allMZ,spmat] = h5watersV2matrix(sp)
%h5watersV2matrix - generate a matrix of intensities following the
%processing in the sp matrix

tic

% Now we need to place into a matrix, which is formed from all values...
all = horzcat(sp{:});
allMZ = unique(all(1,:));

% New matrix
numS = numel(sp);
spmat = sparse(numS,numel(allMZ));

% For each scan whack it in the sparse matrix
for n = 1:numS
    
    tmp = unique(sp{n}','rows');
    
    [~,ia,~] = intersect(allMZ',tmp(:,1));
    spmat(n,ia) = tmp(:,2);
end

toc

end


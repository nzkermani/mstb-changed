function D = distance_mat(L);

%DISTANCE_MAT computes the distance matrix for an n x 2 matrix of
%coordinates,L.

n = length(L);

L = L';

L2 = sum(L.^2,1);

D =sqrt(repmat(L2,n,1) + repmat(L2',1,n) - 2*L'*L);

% Avoid complex-valued zeros

for i = 1:n
    
    D(i,i) = 0;
    
end




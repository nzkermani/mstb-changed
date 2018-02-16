function D = NN_distance_matrix(points);
%%NN_distance_matrix calculates mean of the minimum Euclidean distance between 2D coordinates
%INPUT:points, n*2 coordinates
%OUTPUT:D: distance matrix
%%Author Nazanin Zounemat Kermani, Imperial college London, 2016
n=size(points,1);
points = points';
L2 = sum(points.^2,1);
D =sqrt(repmat(L2,n,1) + repmat(L2',1,n) - 2*points'*points);
%add to diagonal to avoid zero distances of one point to itself
D = D + (max(max(D))* eye(n)) ;
meanNearestNeighbour = mean(min(D));
end


function dist = polyPointDistAndNearestNeighbor(poly,m,N)
%% polyPointDistAndNearestNeighbor simulates the sampling distribution of average 
% minmum nearest-neighbor distance in a fixed polygon. .
% INPUTS:   poly = boundary file of polygon
%           m  = number of points in polygon
%           N  = number of simulations
% OUTPUTS: dist = vector of mean nearest-neighbor distances
%% Author: Nazanin z. Kermani, Imperial college London, 2016
%vector of the average minimum nearest neighbor 
dist = zeros(N,1);
for sim = 1:N
  % Simulate point pattern
  pts = random_location(poly,m);
  % calculate mean nearest neighbors(Euclidean distance)
  dist(sim) =NN_distance_matrix(pts);
end

end

    
    
    
    
    
    
    
    

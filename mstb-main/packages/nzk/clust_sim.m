function OUT = clust_sim(poly,a,m,N,opts,display)
%% CLUST_SIM.M simulates the sampling distribution of average 
% nearest-neighbor distance in a fixed polygon. It can also determine 
% the P-value for a given mean nearest-neighbor distance, if supplied.
% 
% INPUTS: 
%      (i)  poly = boundary file of polygon
%     (ii)    a  = area of polygon
%    (iii)    m  = number of points in polygon
%     (iv)    N  = number of simulations
%      (v)   opts = an (optional) structure with variable inputs:
%                 opts.bins = number of bins in histogram (default = 10)
%                 opts.m_dist = mean nearest-neighbor distance for testing
%                  
% OUTPUTS: OUT = vector of mean nearest-neighbor distances
%
% SCREEN OUTPUT: (i)   mean of simulated mean distances
%                (ii)  histogram of simulated mean distances
%                (iii) P-value of opts.m_dist (if included)
% FUNCTIONS CALLED:  rand_loc.m, polyform.m, distance.m

%% Written by: TONY E. SMITH, 12/31/00
% Modified by Nazanin z. Kermani
if nargin == 5 %parse info if exists 
    
    fields = fieldnames(opts);
    
	nf = length(fields);
	
	if nf > 0 
        
        for i=1:nf
            
            if strcmp(fields{i},'bins')
                bins = opts.bins; 
            elseif strcmp(fields{i},'m_dist')
                d = opts.m_dist;  
            end
        end
        
	else %user input a blank structure
        nf = nf;
	end
    
end %end parsing

%initializations

Z = polyform(poly);

dist = zeros(N,1);

index = 10;

for sim = 1:N
   
   % Simulate point pattern
   
   pts = rand_loc(poly,Z,m);    

	%Construct Nearest-Neighbor Distances

	M = distance_mat(pts);

	area = a;

	mx = max(max(M));

	M = M + mx*eye(m) ; %add to diagonal to avoid zero case

	D = min(M);

	dist(sim) = mean(D);

   
end
OUT = dist;
end

    
    
    
    
    
    
    
    

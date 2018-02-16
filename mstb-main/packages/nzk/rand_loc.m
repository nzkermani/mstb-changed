function LOC = rand_loc(X,Z,N)

%% RAND_LOC distributes random points in a series of polygons

% FUNCTIONS CALLED: rand_num.m, pt_in_poly.m

% INPUTS: X = boundary file in format of 'polyform.m'
%         Z = output of 'polyform.m'
%         N = vector listing number of points for each polygon.

% OUTPUTS: LOC = vector of all points


%out = zeros(1,2); % Place holders for rand_num output
%% Modified by Nazanin z. Kermani for repeatability  


rng(1,'twister');
s = rng;
rng(s);

%number of polygons
k = length(N) ; 

u = ones(k,1) ;

M = N'*u ;

LOC = zeros(M,2) ;

i = 1;

position = 0 ;

while i <= k
   
   if N(i)== 0
      
      i = i  ;      
      
   else
          
   	% Define current polygon
   
   	first = Z(i,1) ;
   	last = Z(i,2) ;
   
   	P = X(first:last,:) ; 
   
   	% Define containing rectangle
   
   	Mn = min(P) ; 
   	Mx = max(P) ;
   
   	x0 = Mn(1) ;
   	y0 = Mn(2) ;
   	x1 = Mx(1) ;
   	y1 = Mx(2) ;
   
   	% Generate points
   
   	j = 1 ;
   
   	while j <= N(i) ;
      
      	kk = 0;
      
      	counter = 1 ;
      
         while (kk == 0)
            
            r = rand(1,1);
                                                      
      		x = x0 + r*(x1 - x0) ;
            
            r = rand(1,1);            
            
            y = y0 + r*(y1 - y0) ;
      
      		if inpolygon(x,y,X(2:end,1),X(2:end,2)) == 1
         
         	LOC(position + j,:) = [x,y] ;
         
         		kk = 1 ;
         
         	end %end if
         
         	% Check to avoid infinite cycles
         
         	counter = counter + 1;
         
         	if (kk == 0)&(counter == 100)
            
            	error('looping error') 
            
         	end	         
         
      	end %end point generation
      
      	j = j + 1 ;
      
   	end 
   
	end % end polygon 
   
   position = position + N(i) ;
   
   i = i + 1 ;
   
end

   
      
         
         
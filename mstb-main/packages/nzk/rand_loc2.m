function LOC = random_location(Poly_,N)

%% random_location generate random points a polygons
%
% INPUT:  Poly_ = boundary file in format of 'polyform.m'
%         N = number of points a polygon.

% OUTPUT: LOC = vector of all points
%% Author Nazanin z. Kermani, Imperial College London 2016

LOC = zeros(N,2) ;
% Define containing rectangle
[xy0 ] = min(Poly_(2:end,:));
[xy1] = max(Poly_(2:end,:));
counter = 1 ;
while counter <= N ;           
            r = rand(1,1);
            x = xy0(1) + r*(xy1(1) - xy0(1)) ;
            r = rand(1,1);            
            y = xy0(2) + r*(xy1(2) - xy0(2)) ;
       		if inpolygon(x(end),y(end),Poly_(2:end,1),Poly_(2:end,2)) == 1
             	LOC(counter ,:) = [x,y] ;
                counter = counter + 1;
            end %end if       
end %end point generation    
end

   
      
         
         
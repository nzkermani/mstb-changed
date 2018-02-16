function  a = area_poly(x,y)

% AREA_POLY calculates the Area of a planar polygon.

%INPUTS: [x,y] = coords of polygon points.

%  Copyright (c) 1995 by Kirill K. Pankratov,
%	kirill@plume.mit.edu.
%	04/20/94, 05/20/95  

%Modified by T.E.Smith Jan.19,2004

% Check Polygon Closure

if (x(1) ~= x(end))
	
    x = [x(:); x(1)];
    
end

if (y(1) ~= y(end))
	
    y = [y(:); y(1)];
    
end

% Calculate contour integral Int -y*dx  (same as Int x*dy).

lx = length(x);

a = -(x(2:lx)-x(1:lx-1))'*(y(1:lx-1)+y(2:lx))/2; 

a = abs(a); %ensures positive area


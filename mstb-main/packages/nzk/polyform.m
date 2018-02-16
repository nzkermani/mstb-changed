function OUT = polyform(X)

% POLYFORM creates first and last rows of boundary file X defining each polygon.

% Written by T.E.Smith, 1/29/99

% INPUT: X = (n:2) matrix describing boundaries of k polygons in the form:
%
%          [1   n1
%           x11 y11 
%           x12 y12
%           .
%           .
%           x11 y11
%           2   n2
%           x21 x22
%           .
%           .
%           .
%           k   nk
%           xk1 yk1
%           .
%           .
%           xk1 yk1]

% OUT: Z = (k:2) matrix listing first and last row of X defining each 
%           polygon.
%        =  [F1  L1  
%           .
%           .
%            Fk  Lk]
%            where, for example, [F1 L1] = [2 n1+1]

dbstop if error;

n = length(X(:,1)) ;


% Count Polygons and check polygon closure

j = 2 ;

kk = 0 ; % Polygon counter

while j < n
   
   k1 = j ;
   
   k2 = X(j-1,2) + (j-1) ; 
   
     
   kk = kk + 1 ;
   
   diff = max(abs(X(k1,:)-X(k2,:)));
   
      
   if (k2<n)&(X(k2+1,1)~=kk+1) %Check polygon numbering
      
      Polygon = kk;
      
      error(['Polygon number ',num2str(kk),' is out of sequence']);
      
   end

   
   if diff > 0 %Check for closure
      
      Polygon = kk
      
      error('POLYGON NOT CLOSED')
      
   end
        
   j = X(j-1,2) + (j+1) ;
   
          
end

% At this point kk = number of polygons

Z = zeros(kk,2);

i = 1 ;

j = 2 ;

while j < n
   
   Z(i,1) = j ;
   
   Z(i,2) = X(j-1,2)+ (j-1) ;   
     
   j = X(j-1,2) + (j+1) ;
   
   i = i + 1 ;
   
end

OUT = Z ;
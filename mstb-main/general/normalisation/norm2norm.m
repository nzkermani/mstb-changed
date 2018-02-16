function [ x] = norm2norm( x )
%norm2norm - normalise a data matrix according to the Euclidean norm of 
% each observation 

numObs = size(x,1);

vect = zeros(numObs,1);

for n = 1:numObs
    
    vect(n,1) = norm(x(n,:),2);
    
end

%tic = nansum(x,2);

x = bsxfun(@rdivide,x,vect);

return

figure; hold on;
scatter(vect,tic,'r');

figure; hold on;
scatter(vect,nansum(x,2),'b');

end


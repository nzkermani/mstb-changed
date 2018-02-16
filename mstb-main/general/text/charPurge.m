function [op] = charPurge(ip,replaceWith)
% charPurge - remove everything that is non-alpha/numeric from a string

% Specify choice of replacement
if nargin == 1
    charRep = '_';
elseif nargin == 2
    charRep = replaceWith;
end

% Find, and replace
lx = isstrprop(ip, 'alphanum');
fx = lx == 0;
ip(fx) = charRep;
op = ip;

end

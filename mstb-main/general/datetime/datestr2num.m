function [ op ] = datestr2num(ip,format)
%datestr2cell - convert a vector of time numbers into a cell format

if nargin == 1
    format = 'yyyy';
end

% Create empty cell
op = zeros(size(ip));

% Demand a certain format
switch format
    case {'dd','mm','yy','yyyy'}
        
    otherwise
        disp('Not numeric');
        return
end

% Now do the rest
for n = 1:size(ip,1)
    
    op(n) = str2double(datestr(ip(n),format));
    
end

end


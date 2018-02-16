function [ op ] = datestr2cell(ip,format)
%datestr2cell - convert a vector of time numbers into a cell format

% Create empty cell
op = cell(size(ip));

% Simple check to make sure format is suitable...
try
    fx = datestr(ip(1),format);
catch
    disp('Unrecognised format - defaulting...');
    format = 'yy/mm/dd';
end

% Now do the rest
for n = 1:size(ip,1)
    
    op{n} = datestr(ip(n),format);
    
end

end


function [op] = statsDateConvert(ip,format)
%statsDateConvert - convert a datestr numeric vector into a cell. Format
%can be specified, or left blank and edited in here on demand...

% Decide default format...
if nargin == 1
    format = 'yyyy';
end
if isempty(format)
    format = 'mm-d';
end

% Convert
op = cell(size(ip));
for n = 1:size(op,1)
    
    op{n,1} = datestr(ip(n),format);
    
end

end


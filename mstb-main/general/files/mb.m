function [ varargout ] = mb(var)
%mb - how large is this variable

w = whos('var');
m = w.bytes / 1e6;

if nargout == 1
    varargout{1} = m;
else
    if m < 1
        disp(['Size: ' sprintf('%0.3f',m*1000) ' Kb']);
    else
        disp(['Size: ' sprintf('%0.3f',m) ' Mb']);
end
    

end


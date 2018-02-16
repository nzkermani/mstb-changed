function [pf,sf] = previousFolder(path)
%previousFolder - return the path of the parent folder

% Decide if it is a string or a cell...
if iscell(path)
    disp('Cell');    
    numP = size(path,1);
    pf = cell(numP,1);
    sf = cell(numP,1);
    
    % Loop through each...
    for n = 1:numP
        [pf{n,1},sf{n,1}] = getInfo(path{n,1});
    end        
    
elseif ischar(path)
    
    % Run the single function here, then return
    [pf,sf] = getInfo(path);
    
else
    disp('CANNOT');
    pf = [];
    sf = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pf,sf] = getInfo(path)

% Add stuff to end...
if ~strcmp(path(end),filesep)
    path = [path filesep];
end

% Find slashes
sl = strfind(path,filesep);

% New path
pf = path(1:sl(end-1));
sf = path(sl(end-1)+1:end-1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
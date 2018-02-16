function [ op ] = statsMZRTRTMZ(mz,rt,method)
%statsMZRTRTMZ - format output in string format according to method 1 or 2,
%where 1 = mzrt and 2 = rtmz

numE = numel(mz);
op = cell(numE,1);

% Which comes first?
switch method    
    case 1
        dt = sortrows([mz rt],1);
    case 2
        dt = sortrows([rt mz],1);
    case 3
        if size(mz,1) > size(mz,2)
            dt = [mz rt];
        else
            dt = [mz' rt'];
        end
end
        
% Loop through
for n = 1:numE
    
    op{n,1} = [sprintf('%0.4f',dt(n,1)) ' | ' sprintf('%0.4f',dt(n,2))];
    
end
    


end


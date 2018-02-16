function [opts] = vararginFcn(argsin)
% Read the arguments and then the data if it wasn't passed

% Define the defaults here
opts.ellipse = true;
opts.ellipseCI = 95;
opts.groups = [];
opts.legend = true;
opts.xLab = 'x-axis';
opts.yLab = 'y-axis';
opts.fsLab = 18;
opts.fsAxes = 16;
opts.markers = {'o','s','d'};
opts.points = true;
opts.cbLab = '';
opts.textLabels = [];
opts.connect = false;
opts.showCB = false;
opts.centroid = false;

% Run through each pair
nArgs = length(argsin);
for i = 1:2:nArgs
    if strcmpi('ellipse',argsin{i})
        tmp = argsin{i+1};        
        if islogical(tmp)
            opts.ellipse = tmp;
        end
        
    elseif strcmpi('confidence',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            if tmp == 90 || tmp == 95 || tmp == 99 || tmp == 99.9
                opts.fP = tmp;
            end
        end
        
    elseif strcmpi('points',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.points = tmp;
        end
        
    elseif strcmpi('centroid',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.centroid = tmp;
        end
        
    elseif strcmpi('groups',argsin{i})
        tmp = argsin{i+1};
        opts.groups = tmp;
        
    elseif strcmpi('legend',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.legend = tmp;
        end
        
    elseif strcmpi('connect',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.connect = tmp;
        end
        
    elseif strcmpi('cb',argsin{i})
        tmp = argsin{i+1};
        if islogical(tmp)
            opts.showCB = tmp;
        end
        
    elseif strcmpi('xlabel',argsin{i})
        tmp = argsin{i+1};
        if ischar(tmp)
            opts.xLab = tmp;
        end
    elseif strcmpi('ylabel',argsin{i})
        tmp = argsin{i+1};
        if ischar(tmp)
            opts.yLab = tmp;
        end
        
    elseif strcmpi('cblabel',argsin{i})
        tmp = argsin{i+1};
        if ischar(tmp)
            opts.cbLab = tmp;
        end
        
    elseif strcmpi('sizelabel',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            opts.fsLab = tmp;
        end
        
    elseif strcmpi('sizeaxes',argsin{i})
        tmp = argsin{i+1};
        if isnumeric(tmp)
            opts.fsAxes = tmp;
        end
        
    elseif strcmpi('confidence',argsin{i})
        tmp = argsin{i+1};
        opts.markers = tmp;       
    
    elseif strcmpi('text',argsin{i})
        tmp = argsin{i+1};
        opts.textLabels = tmp;       


    
    end
end

end

function [locn] = statsLocations(layout,type)
%desiLocations - these are the locations for the various axes and labels
%and things in the desi window...

% We have a few things to draw, many of which can just have their
% visibility changed as we switch the views...
switch lower(type)
    case 'ms'
        [locn] = defMS;
    case 'lcms'
        [locn] = defLCMS;
end

% Which of the layouts (i.e. table or axes) to show?
switch layout
    
    case 'axes'
        
        locn.tab.vis = 'off';        
        
    case 'table'
        
        locn.scatter.vis = 'off';
        locn.conf.vis = 'off';
        locn.load.vis = 'off';
        if strcmp(type,'lcms')
            locn.spec.vis = 'off';
        else
            locn.spec.vis = 'on';
        end
        
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [locn] = defMS
% Default positions for MS

% Observation table
locn.tab.pos = [0.01 0.31 0.98 0.68];
locn.tab.vis = 'on';

% Spectral plot
locn.spec.pos = [0.01 0.05 0.98 0.13];
locn.spec.vis = 'on';

% Scatter axes
locn.scatter.pos = [0.01 0.41 0.48 0.58];
locn.scatter.vis = 'on';

% Confusion axes / box plot etc...
locn.conf.pos = [0.51 0.41 0.48 0.58];
locn.conf.vis = 'on';

% Loading axes
locn.load.pos = [0.01 0.21 0.98 0.18];
locn.load.vis = 'on';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [locn] = defLCMS
% Default positions for LC-MS

% Observation table
locn.tab.pos = [0.01 0.31 0.98 0.68];
locn.tab.vis = 'on';

% Spectral plot
locn.spec.pos = [0.51 0.05 0.48 0.42];
locn.spec.vis = 'on';

% Scatter axes
locn.scatter.pos = [0.01 0.51 0.48 0.48];
locn.scatter.vis = 'on';

% Confusion axes / box plot etc...
locn.conf.pos = [0.51 0.51 0.48 0.48];
locn.conf.vis = 'on';

% Loading axes
locn.load.pos = [0.01 0.05 0.48 0.42];
locn.load.vis = 'on';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

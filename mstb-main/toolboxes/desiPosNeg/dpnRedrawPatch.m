function dpnRedrawPatch(~,~,fig,type)
%dpnRedrawPatch - re draw the designated annotations, i.e. this function
%needs to be called for each of the three axes as required

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Check that there are some annotations
if ~isfield(dpn,'anno')
    return
end
numA = size(dpn.anno,1);


% Which annotations and where and which axes
switch type
    case 'opt'
        
        % Define the axes
        axes(dpn.fig.ax.opt(1));
        hold on;
        
        % Loop through...
        for n = 1:numA    
            x = dpn.anno{n,6};
            y = dpn.anno{n,7};
            x = [x(1) x(2) x(2) x(1)];
            y = [y(1) y(1) y(2) y(2)];
            if x(1) == x(2) && y(1) == y(3)
                x = x(1);
                y = y(1);
            end
            dpn.anno{n,1}(1) = redrawPatch(x,y,dpn.anno{n,3});
        end
        
    case 'ms1'
        % MS1
        axes(dpn.fig.ax.ms1(1));
        hold on;
        for n = 1:numA        
            dpn.anno{n,1}(2) = redrawPatch(dpn.anno{n,8},...
                dpn.anno{n,9},dpn.anno{n,3});
        end

    case 'ms2'
        % MS2
        axes(dpn.fig.ax.ms2(1));
        hold on;
        for n = 1:numA        
            dpn.anno{n,1}(3) = redrawPatch(dpn.anno{n,8},...
                dpn.anno{n,9},dpn.anno{n,3});
        end

    otherwise
        error('No such things');
end

% Update as we need to save the handles to the new patches...
guidata(fig.fig,dpn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h] = redrawPatch(x,y,col)
% Draw a patch in the CURRENT axes


if numel(x) == 1 && numel(y) == 1
    % Then this is a single scatter point of annotation
    h = scatter(x(1),y(1),200,col,'o','filled',...
        'MarkerEdgeColor','w');
    
else
    
    % Check that we have the annotations in the correct order for optical
    % and MS annotations
    x = [min(x) max(x) max(x) min(x)];
    y = [min(y) min(y) max(y) max(y)];

    % Then multiple points, so we patch them together
    x = x + [-0.5 0.5 0.5 -0.5];
    y = y + [-0.5 -0.5 0.5 0.5];
    
    h = patch(x,y,...
        col,...
        'EdgeColor',col,...
        'FaceColor',col,...
        'FaceAlpha',0.4,...
        'LineWidth',3);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

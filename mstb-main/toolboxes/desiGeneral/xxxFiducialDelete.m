function xxxFiducialDelete(src,event,fig,man)
%xxxFiducialDelete - remove the selected markers

% Guidata
dpn = guidata(fig.fig);

% Get the table...
if isfield(dpn.opt,'fiducial')
    if isempty(dpn.opt.fiducial)
        return
    else
        tab = get(man.tab,'Data');
    end
else
    return
end

% Find out which ones have been clicked...
fx = cell2mat(tab(:,1)) == true;
if sum(fx) == 0
    return;
end

% Now we must delete the scatter points from the axes
delete([dpn.opt.fiducial{fx,5}]);
delete([dpn.opt.fiducial{fx,6}]);

% Now trim out the deleted ones...
dpn.opt.fiducial = dpn.opt.fiducial(~fx,:);
set(man.tab,'Data',dpn.opt.fiducial(:,1:2));

% Update metadata
guidata(fig.fig,dpn);


end


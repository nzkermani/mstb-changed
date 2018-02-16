function dpnInterpolateCallback(~,~,fig)
%dpnInterpolateCallback - function that calls the interpolation method,
%which is useful for performing interpolation using different settings

% Get the guidata and return
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Ask the user for the method to be performed...
allM = {'linear','cubic','spline','nearest','bilinear','bicubic'};
[intMethod,ok] = listdlg('ListString',allM,...
    'SelectionMode','single',...
    'Name','Interpolation Method');

% See if anything was selected
if ~ok
    return
end

% This is the interpolation method, which we will perform for both of the
% ion modes in the same way
intMethod = allM{intMethod};

% For each of the modes, we need to determine the vector that says which
% pixels are original and which are interpolated. Note that the variable we
% use is isInterp, to which we want the opposite values.
i1 = dpn.d1.isInterp(:,1) == 0;
i2 = dpn.d2.isInterp(:,1) == 0;

% Interpolation D1
[dpn.d1.sp,...
    dpn.opts.interp,...
    dpn.d1.isInterp] = dpnPerformInterpolation(dpn.d1.sp,intMethod,i1);

% Interpolation D2
[dpn.d2.sp,...
    dpn.opts.interp,...
    dpn.d2.isInterp] = dpnPerformInterpolation(dpn.d2.sp,intMethod,i2);

% Save to the guidata
guidata(fig.fig,dpn);

% Finally update the ion images...
dpnIonImage([],[],fig.ax.ms1,nansum(dpn.d1.sp,3));
dpnIonImage([],[],fig.ax.ms2,nansum(dpn.d2.sp,3));

end


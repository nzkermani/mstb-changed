function dpnComparisonCallback(~,~,fig)
%desiComparisonCallback - perform all interpolations, match the peaks, and
%then output the results. Be careful not to save any of the data, as
%bilinear and bicubic implementations change the original data

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Which methods?
methods = {'linear','cubic','spline','nearest','bilinear','bicubic'};
numM = numel(methods);

% For each of the modes, we need to determine the vector that says which
% pixels are original and which are interpolated. Note that the variable we
% use is isInterp, to which we want the opposite values.
i1 = dpn.d1.isInterp(:,1) == 0;
i2 = dpn.d2.isInterp(:,1) == 0;

% Loop through the various interpolations
for n = 1:numM
        
    % Clear new, then add mz
    clear new
    new.d1.mz = dpn.d1.mz;
    new.d2.mz = dpn.d2.mz;    
    
    % Interpolate d1
    [new.d1.sp,...
        new.d1.isInterp,~] = dpnPerformInterpolation(dpn.d1.sp,methods{n},i1);
    
    % Interpolate d2
    [new.d2.sp,...
        new.d2.isInterp,~] = dpnPerformInterpolation(dpn.d2.sp,methods{n},i2);
    
    % Run the comparison function
    [op] = dpnPaperComparison(new,true);
    
    % Now we need to extract the various bits and pieces
    crr.(methods{n}) = op.crr;   
    
end

% Return crr to the workspace
assignin('base','crr',crr);











return

% What is the interpolation method?
im = dpn.opts.interp;
disp(['INTERP method = ' im]);

% Run the other function that does peak matching etc
[op] = dpnPaperComparison(dpn,true);

% How about just exporting the correlations to the workspace, using the
% name of the correlation... First get the names of the variables to see if
% there is already the one we wnat
allV = evalin('base','who');
fx = strcmp(allV,'crr');
if sum(fx) == 1
    crr = evalin('base','crr');
end

% Add to the structure
crr.(im) = op.crr;
assignin('base','crr',crr);

disp('-');

end


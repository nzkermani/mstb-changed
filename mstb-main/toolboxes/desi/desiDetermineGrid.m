function desiDetermineGrid(~,~,fig)
%desiDetermineGrid - calculate the grid spacings and what not such that we
%can draw a grid and put the annotations in the right place when it comes
%to that part of the analysis.  Ideally we run this function automatically
%after the co-registration procedure, saving the results to the guidata
%structure

% Get guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Not even an optical image?
if ~isfield(dpn,'opt')
    disp('Upload a proper optical image first');
    return
end

% Ensure that we have done co-registration too
if ~isfield(dpn.opt,'coreg')
    disp('Need to co-register with an optical image first');
    return
end

% Grid determination function
[dpn.fig.grid.opx,dpn.fig.grid.opy,...
    dpn.fig.grid.msx,dpn.fig.grid.msy] = xxxGridCalculate(dpn.opt.coreg,dpn.d1.img);

% Save the grids in the guidata structure
guidata(dpn.fig.fig,dpn);

return

% % What is the size of the opt image?
% szOp = size(dpn.opt.coreg);
% 
% % What is the size of the ms image
% szMS = size(dpn.d1.img);
% 
% % Determine the ratio between the two
% ratio = szOp(1:2) ./ szMS(1:2);
% 
% % Linear interpolation across the optical image. Also adjust the gridlines
% % to intersect neighbouring pixels, rather than intersecting within a pixel
% dpn.fig.grid.opx = (1:ratio(2):szOp(2)) - (ratio(2) / 2);
% dpn.fig.grid.opy = (1:ratio(1):szOp(1)) - (ratio(1) / 2);
% 
% % Determine the MS grid
% dpn.fig.grid.msx = (1:szMS(2))-0.5;
% dpn.fig.grid.msy = (1:szMS(1))-0.5;

end


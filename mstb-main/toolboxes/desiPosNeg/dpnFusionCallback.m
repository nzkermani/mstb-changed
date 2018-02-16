function dpnFusionCallback(~,~,fig)
%dpnFusionCallback - prep the data for the fusion callback window

% Get the guidata and decide what to do
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% What inputs do we need to provide? A structure containing the relevant
% information, such as mz|sp|size|meta. The data should also be in a 2D
% format, so images need to be reshaped...
data = struct('name',[],'x',[],'y',[],'imSize',[],'meta',[],'tobg',[]);
for n = 1:2
    
    % One dataset or the other
    if n == 1
        dx = dpn.d1;
        data(n).name = 'Positive';
        if isfield(dx,'tobg')
            tobg = dx.tobg;
        end
    else
        dx = dpn.d2;
        data(n).name = 'Negative';
        if isfield(dx,'tobg')
            tobg = dx.tobg;
        end
    end
    
    % Size?
    sz = size(dx.sp);
    data(n).imSize = sz;
    
    % Reshape
    data(n).y = reshape(dx.sp,[sz(1)*sz(2) sz(3)]);
    data(n).x = dx.mz;
    
    % Sundries... to be added here when they may be necessary
    
end

% Need to get the annotations from the data
[anno.mask2,...
    anno.histID,...
    anno.isInt1,...
    anno.isInt2,...
    anno.pixID] = dpnAnnotationExtract(dpn);

% Close out the tobg image to make it intact
mOpts.imopen = 2;
mOpts.imclose = 8;
mOpts.bigoper = [];
tobg = doMorph(tobg,mOpts);

data(1).tobg = tobg;
data(2).tobg = tobg;

% Run the function...
dataFusion(data,anno,dpn);

end


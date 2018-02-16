function statsWorkspace(~,~,fig)
%statsWorkspace - get variables directly from the workspace...

% Clear all existing axes
statsAxesClearAll([],[],fig);

% Get the names of the variables...
varNames = evalin('base','who');
if isempty(varNames)
    return
end

% Present names as a listbox
[sel,~] = listdlg('PromptString','Select variable',...
    'SelectionMode','single',...
    'ListString',varNames);

% Try to import that which was selected
tmp = evalin('base',varNames{sel});

% Let's consider it as either an MS structure or an LCMS structure
if isfield(tmp,'v') && isfield(tmp,'x') && isfield(tmp,'y')
    datatype = 'lcms';
elseif isfield(tmp,'mz') && isfield(tmp,'sp') && isfield(tmp,'meta')
    datatype = 'ms';
else
    datatype = 'none';
end

% Now import the dataset
switch datatype
    
    case 'ms'
        
        if size(tmp.mz,2) > size(tmp.mz,1)
            sts.raw.var.mz = tmp.mz';
        else
            sts.raw.var.mz = tmp.mz;
        end
        sts.raw.sp = tmp.sp;
        sts.raw.meta = tmp.meta;
                
    case 'lcms'
        
        sts.raw.var = tmp.v;
        sts.raw.sp = tmp.x;
        sts.raw.meta = tmp.y.meta;
        
        % There may be other meta variables to be introduced
        
        % Disable some functions that don't yet work for LC-MS data
        fig.tb.doRatio.Enable = 'off';
        fig.tb.doExternal.Enable = 'off';
        
    otherwise
        
        disp('Cannot import data structure');
        return
end

% NEW: add a numeric vector from 1:n to track obs from raw->proc in meta.
% This is going to have to be a hidden thing as it is mostly useless.
vec = 1:size(sts.raw.sp,1);
sts.raw.obsID = vec';

% Save into a new sts guidata structure
sts.proc = sts.raw;
sts.fig = fig;
sts.datatype = datatype;
sts.res = [];

% What about the path?
if isfield(tmp,'path')
    sts.path = tmp.path;
else
    sts.path = [];
end

% Save as guidata
guidata(fig.fig,sts);

% Change the layout to match what it should be...
statsLayoutChange([],[],fig,datatype);

% Now create a function that plots the data in the spectral window
statsPlotSpectra([],[],fig,'raw');

% Populate the table...
statsTablePopulate([],[],fig);

end


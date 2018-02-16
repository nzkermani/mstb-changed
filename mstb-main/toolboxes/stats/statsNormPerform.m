function statsNormPerform(~,~,fig,window)
%statsNormPerform - do something to the data...

% Get the guidata
sts = guidata(fig.fig);
if isempty(sts)
    return
end

% Get the observations from the table which are still ticked. These are the
% only ones that we care for...
td = get(fig.ax.table,'Data');
sel = cell2mat(td(:,1));

% Get the variable selection function...
varSel = window.varSelect.String{window.varSelect.Value};
switch varSel
    
    case 'Developmental'
    
        % Run a function here to check it all out...
        varMask = statsVariableFiltering(sts);

    case 'Frequent'
        
    case 'Very frequent'
        freq = sum(sts.raw.sp > 0,1) / size(sts.raw.sp,1);
        varMask = freq >= 0.05;
        
    otherwise
        % Do nothing, i.e. all variables
        varMask = sts.raw.var.mz' >= 0;
end

% Trim the metadata as well
fn = fieldnames(sts.raw.meta);
for n = 1:numel(fn)
    tmp = sts.raw.meta.(fn{n});
    meta2.(fn{n}) = tmp(sel,:);
end

% Run the function to change the data
if strcmp(sts.datatype,'ms')
    [sts.proc.var.mz,sts.proc.sp,sts.proc.opts,~] = xxxNormTranScal(window,...
        sts.raw.var.mz(varMask),...
        sts.raw.sp(sel,varMask),...
        meta2);

elseif strcmp(sts.datatype,'lcms')
    [tmp,sts.proc.sp,sts.proc.opts,mask] = xxxNormTranScal(window,...
        [sts.raw.var.mz(varMask) sts.raw.var.rt(varMask)],...
        sts.raw.sp(sel,varMask));
    sts.proc.var.mz = tmp(:,1);
    sts.proc.var.rt = tmp(:,2);
    sts.proc.var.rtmz = sts.raw.var.rtmz(mask);
    sts.proc.var.all = sts.raw.var.all(mask,:);
    
end

% Check for failure
if isempty(sts.proc.opts)
    return
end

% Trim the metadata as well
meta = sts.raw.meta;
fn = fieldnames(meta);
for n = 1:numel(fn)
    tmp = meta.(fn{n});
    meta.(fn{n}) = tmp(sel,:);
end
sts.proc.meta = meta;
sts.proc.obsID = sts.raw.obsID(sel,1);

% As we have changed the data, we need to delete any existing results that
% have been stored already
sts.res = [];

% Update the guidata
guidata(fig.fig,sts);

% Update the plot at the bottom...
statsPlotSpectra([],[],fig,'proc');

end


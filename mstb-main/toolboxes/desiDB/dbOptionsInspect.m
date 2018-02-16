function dbOptionsInspect(src,event,data,fig,pan)
%dbOptionsInspect - determine the processing options and prepare the data
%for analysis

% Get the various options from the window
tmp = get(pan.pixel,'String');
val = get(pan.pixel,'Value');
opts.pixels = tmp{val};

% Normaliation
tmp = get(pan.norm,'String');
val = get(pan.norm,'Value');
opts.norm = tmp{val};

% m/z range
tmp = get(pan.mzRange,'String');
numR = size(tmp,1);
opts.mzRange = zeros(numR,2);
for n = 1:numR
    
    val = str2num(tmp(n,:)); %#ok<ST2NM>
    if numel(val) == 2
        opts.mzRange(n,:) = val;
    elseif numel(val) == 1
        opts.mzRange(n,:) = [val val+1];
    end    
end
fx = sum(opts.mzRange,2) > 0;
opts.mzRange = opts.mzRange(fx,:);
set(pan.mzRange,'String',num2str(opts.mzRange));
    
% Now the stuff about the table and re-annotations...
td = get(pan.table,'Data');
grps = get(pan.table,'ColumnName');
orig = get(pan.table,'RowName');
new = repmat({'Delete'},[numel(orig) 1]);

% Decide on the new labelling
for n = 1:size(td,1)
    
    % Ensure that we have only one match
    ticks = cell2mat(td(n,:));
    if sum(ticks) == 1        
        new{n,1} = grps{ticks'};
    else
        % Just delete it...
        td(n,:) = {false};
    end
    
end

% Update the table in case of issues, i.e. double ticking
set(pan.table,'Data',td);

% Save the group names into the options structure
opts.names = [orig new];

% Now we need to launch another function that does the processing, and then
% leads on to the visualisation
dbFunctionInspect([],[],data,opts,fig,pan);

end


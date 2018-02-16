function statsAddMoreMeta(~,~,fig,~)
%statsAddMoreMeta - import stuff from a csv/xlsx file

% Get the data, and hopefully a folder...
sts = guidata(fig.fig);

% If there is a path, then use it
if ~isempty(sts.path)
    defP = sts.path;
else
    defP = [pwd filesep];
end

% Ask the user for a file...
[a,b,~] = uigetfile({'*.csv; *.xlsx; *.xls'},'Load a File',defP);
if isnumeric(a)
    return
end
fileExt = fileExtension(a);

% Treat differently if xlsx or csv?
switch fileExt
    case 'csv'        
        
        % Read in as a table
        dt = readtable([b a],'Delimiter',',');
        
        % Now need to convert it slightly...
        heads = dt.Properties.VariableNames;
        
        
        isNum = false(size(heads));
        numH = numel(heads);
        meta = cell(size(dt,1),numH);
        for n = 1:numH
            try
                meta(:,n) = dt.(heads{n});
            catch
                tmp = dt.(heads{n});
                for r = 1:size(dt,1)
                    meta{r,n} = tmp(r);
                end
                isNum(n) = true;
            end
                
        end
        
    case {'xls','xlsx'}
        [~,~,dt] = xlsread([b a]);
        
        % Convert...
        heads = dt(1,:);
        
        % Logical test to ensure that we don't keep empty heads...
        emptyHeads = ~cellfun(@isempty,...
            cellfun(@true,...
            cellfun(@isnan,heads',...
            'UniformOutput',false),...
            'UniformOutput',false));
        heads = heads(~emptyHeads);
        dt = dt(:,~emptyHeads);

        
        isNum = false(size(heads));
        numH = numel(heads);
        meta = cell(size(dt,1)-1,numH);
        for n = 1:numH
            meta(:,n) = dt(2:end,n);
            
            % Convert to suitable name...
            sp = ~isstrprop(heads{n},'alphanum');
            heads{n}(sp) = '_';
            
            % Work out if it is numeric?
            tt = cellfun(@isnumeric,meta(:,n));
            if sum(tt) == numel(tt)
                isNum(n) = true;
            
            elseif sum(tt) > 0
                % Try to convert the numeric parts into text parts
                ftt = find(tt);
                for r = 1:numel(ftt)
                    meta{ftt(r),n} = num2str(meta{ftt(r),n});
                end
                
            end
            
        end       
        
end

% Remove trailing and leading ' marks from the meta(:,1)
fx = strfind(meta(:,1),'''');
for n = 1:size(fx,1)
    if ~isempty(fx{n,1})
        if numel(fx{n,1}) >= 2
            meta{n,1} = meta{n,1}(2:end-1);
        end
    end
end
    

% Now add the stuff to the table by matching file names

% Find unique entries in sts.raw.meta.fileID
[unq,~,ind] = unique(sts.raw.meta.fileID);
numE = numel(unq);


% Do separately for each header
for h = 2:numH
    
    % Storage
    tmp = cell(numel(ind),1);
    
    % Loop through each entry
    for n = 1:numE
        
        % Find matches to unq{n}
        fx = strcmpi(meta(:,1),[unq{n} '.mat']);
        if sum(fx) == 0
            fx = strcmpi(meta(:,1),unq{n});
        end
        
        % Add bits to the right place
        fy = ind == n;
        
        % What if there is no match?
        if sum(fx) > 0
            tmp2 = meta(fx,h);
            tmp(fy,1) = tmp2(1);
        else
            if ~isNum(h)
                tmp(fy,1) = repmat({'Unknown'},[sum(fy) 1]);
            else
                tmp(fy,1) = repmat({NaN},[sum(fy) 1]);
            end
        end
        
    end
    
    % Now that it is finished, we save it to the structure
    if ~isNum(h)
        sts.raw.meta.(heads{h}) = tmp;
    else
        sts.raw.meta.(heads{h}) = cell2mat(tmp);
    end
    
end

% Now let's update the thing and be done with it...
guidata(fig.fig,sts);

% Then we need to update the table?
statsTablePopulate([],[],fig);
  
end


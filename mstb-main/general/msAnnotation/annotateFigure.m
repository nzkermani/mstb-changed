function annotateFigure(list,annoInfo,annoStr,type)
%annotateFigure - plot something useful to show which part is
%interesting...

% Determine if 'type' is a category or subcategory
sz = size(annoStr,2);
category = cell(sz,1);
for n = 1:sz
    category{n,1} = annoStr(n).category;
end

% Determine if if was one of these categories which was picked
fx = strcmpi(category,type);
if sum(fx) == 1

    % Plot the appropriate bar chart (perhaps other figure sometime?)
    data = annoStr(fx).table;
    rowName = annoStr(fx).values;
    colName = annoStr(fx).histID;
    
    figure('Units','normalized','Position',[0.25 0.5 0.5 0.4]);
    h = bar(data);
    set(gca,'XTickLabels',rowName,'XTick',1:numel(rowName),'FontSize',14);
    legend(h,colName);
    ylabel('Number of significant features','FontSize',18,'FontWeight','bold');

    return
end

% If neither of these then we need to determine which subgroup it was. So
% we need a list of all the categorical subgroups...
subgrp = vertcat(annoStr.values);

% Extract the numeric bits
isnu = cellfun(@isnumeric,subgrp);
subNum = NaN(numel(isnu),1);
subNum(isnu,1) = cell2mat(subgrp(isnu));

% Extract the text bits
subTxt = cell(numel(isnu),1);
subTxt(~isnu,1) = subgrp(~isnu);

% Decide which is the match
if isnumeric(type)
    fx = subNum == type;
else
    fx = strcmpi(subTxt,type);
end

% Ensure that we have a match
if sum(fx) == 0
    disp('No match');
    return
elseif sum(fx) > 1
    disp('Multiple matches?');
    return
end

% So now that we have a match we look at finding those entries in the
% larger list... But we need to determine what it is...
qry = subgrp(fx);
if isnumeric(type)
    qry = cell2mat(qry);
end

% Find it in the structure
idx = false(sz,1);
if isnumeric(type)
    for n = 2:3
        tmp = cell2mat(annoStr(n).values) == qry;
        if sum(tmp) == 1
            idx(n,1) = true;
        end
    end
else
    for n = [1 4 5]
        tmp = strcmpi(annoStr(n).values,qry);
        if sum(tmp) == 1
            idx(n,1) = true;
        end
    end
end

% Find the indices of the non empty rows
fx = cellfun(@isempty,annoInfo(:,idx));
gx = find(~fx);

% Now match the
if isnumeric(qry)
    fy = cell2mat(annoInfo(~fx,idx)) == qry;
else
    fy = strcmpi(annoInfo(~fx,idx),qry);
end
  
% Return to a full sized vector
fsv = zeros(numel(fx),1);
fsv(gx,1) = fy; %#ok<FNDSB>

% Determine the specific variables and their histologically significant
% group
table = annoInfo(fsv == 1,:);
subl = list(fsv == 1,:);

% Unique histological groups for these variables
[unqHist,~,unqInd] = unique(subl(:,5));
numG = numel(unqHist);

% Make a table for each of the columns
bpd = struct('rowL',[],'data',[]);
empty = zeros(size(table,2),1);
for n = 1:size(table,2)

    % Determine unique values in this column
    ex = ~cellfun(@isempty,table(:,n));
    tableData = table(ex,n);
    if n == 2 || n == 3
        [unq,~,~] = unique(cell2mat(tableData));
    else
        [unq,~,~] = unique(tableData);
    end
    numL = numel(unq);
    count = zeros(numL,numG);

    % Do also for each histological group
    for g = 1:numG        
        gidx = unqInd(ex) == g;
        
        % For each of the subbits
        for r = 1:numL
            
            if n == 2 || n == 3
                fx = cell2mat(tableData(gidx,1)) == unq(r);
            else
                fx = strcmpi(tableData(gidx,1),unq{r});
            end
            
            count(r,g) = sum(fx);
            
        end            
        
    end
    
    % Save the data for box plots or similar
    bpd(n).rowL = unq;
    bpd(n).data = count;
    if sum(ex) == 0
        empty(n,1) = 1;
    end
    
end

% Now make a nice set of box plots...
bpd(empty == 1) = [];
numE = size(bpd,2)

figure;
for n = 1:numE
    
    subplot(1,numE,n);
    
    if size(bpd(n).data,1) == 1
        tmp = [bpd(n).data; zeros(1,numel(bpd(n).data))];
        
    else
        tmp = bpd(n).data;
    end
    
    h = bar(tmp);
    set(gca,'XTickLabels',bpd(n).rowL,'XTick',1:numel(bpd(n).rowL),'FontSize',14);
    legend(h,unqHist);
    ylabel('Number of significant features','FontSize',18,'FontWeight','bold');

end  
    
        

end


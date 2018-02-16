function [trim] = annotateLipidClass(mz,sp,grp,ass,db)
%annotateOP2 - revised function, hopefully simpler...

% Threshold for significance
qThresh = 0.05;

% Grouping variables
[unq,~,unqInd] = unique(grp);

% Determine the p/q values
[allStat] = determinePQvalues(sp,grp);

% Somewhere to store a list of variables
numV = numel(mz);
anno = cell(numV,7);

% Determine the number of annotations made for each file
totAnno = cellfun(@size,ass.match,'UniformOutput',false);
totAnno = cellfun(@max,totAnno);

% Loop through the variables
numV = size(sp,2);
numAdd = size(ass.match,2);
for n = 1:numV
                 
    % Just combine the indices
    tmp = ass.match(n,:)';
    match = cell2mat(tmp);
    
    % What are the non-empty ones?
    if isempty(match)
        anno(n,:) = {n mz(n) 1 1 'N/A' [] []};
        continue;
    end
    
    % Determine the information about these annotations
    [name,class] = genAnnoInfo(match,db);
   
    % Format the p/q/diff
    p = allStat(n,1);
    q = allStat(n,3);
    if isnan(p)
        p = 1;
    end
    if isnan(q)
        q = 1;
    end
    dfg = unq{allStat(n,2)};
    
    % Place the information into the anno cell matrix
    anno(n,:) = {n mz(n) p q dfg name class};
    
    % Should we store the adduct?
end

% So with an annotated list of significant variables, we should consider
% splitting the variables based on differences. First, let's trim out the
% variables which aren't significant
fx = cell2mat(anno(:,4)) < qThresh;
trim = anno(fx,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [name,class] = genAnnoInfo(match,db)

numM = numel(match);
name = cell(1,numM);
class = cell(1,numM);

for n = 1:numM
    
    % Annotation name
    if isfield(db,'Name')
        tmpName = db.Name{match(n)};
    elseif isfield(db,'Formula')
        tmpName = db.Formula{match(n)};
    else
        tmpName = '?';
    end
    name{1,n} = tmpName;
    
    % Class
    class{1,n} = [db.Class1{match(n)} ',' db.Class2{match(n)} ',' db.Class3{match(n)}];
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avgs,allFC] = genFCs(grp,idx,sp)
% Determine mean/median and then fold changes between groups

numG = numel(grp);
avgs = zeros(numG,2);

% Determine mean / median
for n = 1:numG    
    fx = idx == n;
    avgs(n,:) = [nanmean(sp(fx)) nanmedian(sp(fx))];    
end

% How about fold changes? Need to make it calculate FCs...
numFC = numG * (numG - 1) / 2;
allFC = zeros(1,numFC);
i = 0;
for n = 1:numG    
    for r = n+1:numG
        i = i + 1;
        allFC(i) = avgs(n,2) / avgs(r,2);
    end
end
allFC = log2(allFC);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allStat] = determinePQvalues(sp,grp)
% In order to do correction of the p values into q values, we have to
% calculate the p-values all together initially...

% Store the data here
numV = size(sp,2);
allStat = zeros(numV,3);

% Unique groups
[unq,~,unqInd] = unique(grp);
if numel(unq) == 2
    iX = unqInd == 1;
    iY = unqInd == 2;
end

warning off all

% Loop
for n = 1:numV
   
    % Get the statistics for this variable, either MWU or ANOVA (n > 2)
    if numel(unq) == 2

        % Do mannU
        [allStat(n,1)] = ranksum(sp(iX,n),sp(iY,n));
        
        % Determine the group in which it is biggest
        if median(sp(iX,n)) > median(sp(iY,n))
            allStat(n,2) = 1;
        else
            allStat(n,2) = 2;
        end
        
    else
        % Do ANOVA
        %[a,b] = anovaPH(op.cmz(n),op.XPeaks(:,n),op.histID);
        [a,~,c] = anova1(sp(:,n),grp,'off');
        
        allStat(n,1) = a;
        
        % Find largest group from c.means
        [~,idx] = max(c.means);
        txtGrp = c.gnames{idx};
        idx = strcmp(unq,txtGrp);
        allStat(n,2) = find(idx);
        
    end
    
end

warning on all

% Here we do the q-value correction...
allStat(:,3) = getBHYqVls(allStat(:,1)',0.05);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

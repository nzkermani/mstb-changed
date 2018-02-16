function annotateOP2(mz,sp,grp,ass,db,txtP,qThresh)
%annotateOP2 - revised function, hopefully simpler...

% Here is where we would have all of the pre-amble and file creation. This
% can come once we have the full function working...
[fid,wb] = preamble(txtP,grp,ass);

% Convert qThresh to a number
if ~isnumeric(qThresh)
    switch qThresh
        case 'All'
            qThresh = 1;
        otherwise
            qThresh = str2double(qThresh);
    end
end

% Grouping variables
[unq,~,unqInd] = unique(grp);

% Determine the p/q values
[allStat] = determinePQvalues(sp,grp,wb,qThresh);
%allStat(isnan(allStat)) = 1;

% Loop through the variables
numV = size(sp,2);
numA = size(ass.match,2);
for n = 1:numV
    
    % Update the waitbar
    if mod(n,100) == 0
        waitbar(n/numV,wb,[int2str(n) '/' int2str(numV)]);
    end

    % Determine if the q value is below the threshold set by the user
    if allStat(n,3) > qThresh 
        continue;
    end    
    
    % Print the information about the actual variable, e.g. m/z
    fprintf(fid,'%0.4f\t',mz(n));
                    
    % What are the matches...
    for r = 1:numA
        
        % Match index stored here
        match = ass.match{n,r};
    
        % Skip following part if no match
        if isempty(match)
            %annoMZ = '';
            name = '';
            class = '';
            mzDB = '';
            ppmDB = '';
        else
            
            % Predetermine the m/z values of the annotations
            annoMZ = mass2mz(db.Mass(match),ass.adduct);
        
            % Prepare all the parts of the annotation, such as m/z, text
            % strings and the stuff... Note that there will be one of these
            % for each of the adducts
            [name,class,mzDB,ppmDB] = genAnnoInfo(match,db,annoMZ(:,r),mz(n));
        end
        
        % Print to the file...
        fprintf(fid,'%s\t%s\t%s\t%s\t',name,mzDB,ppmDB,class);       

    end
    
    % Print p/q values
    if isnan(allStat(n,1))
        fprintf(fid,'\t\t\t');
    else
        fprintf(fid,'%0.2E\t%0.2E\t%s\t',allStat(n,1),allStat(n,3),unq{allStat(n,2)});
    end
    
    % Determine mean/median / fold changes
    [avgs,allFC] = genFCs(unq,unqInd,sp(:,n));
    for r = 1:size(avgs,1)
        fprintf(fid,'%0.4f\t%0.4f\t',avgs(r,1),avgs(r,2));
    end
    for r = 1:numel(allFC)
        fprintf(fid,'%0.2f\t',allFC(r));
    end
        
    fprintf(fid,'\n');

end

fclose(fid);
delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fid,wb] = preamble(txtP,grp,ass)

% Define the filename
if isempty(txtP)
    txtP = [datestr(now,'yymmdd-HHMMSS') '.txt'];
end
fid = fopen(txtP,'w');

% Create a waitbar
wb = waitbar(0,'Initiating','Name','Exporting Annotations');

% m/z
fprintf(fid,'m/z\t');

% Information about each adduct
numA = numel(ass.adduct.name);
for n = 1:numA
    fprintf(fid,'%s\t',[ass.adduct.name{n} ' Name']);
    fprintf(fid,'%s\t',[ass.adduct.name{n} ' m/z']);
    fprintf(fid,'%s\t',[ass.adduct.name{n} ' ±ppm']);
    fprintf(fid,'%s\t',[ass.adduct.name{n} ' Class']);
end

% Groups
[unq,~,~] = unique(grp);
if numel(unq) == 2
    stat1 = 'MannU pVal';
    stat2 = 'MannU qVal';
else
    stat1 = 'ANOVA pVal';
    stat2 = 'ANOVA qVal';
end

% p/q values
fprintf(fid,'%s\t%s\t%s',stat1,stat2,'Diff');
   
% Add in the headings for the specific groups...
for n = 1:numel(unq)
    fprintf(fid,'\t%s',[unq{n,1} ' (Mean)']);
    fprintf(fid,'\t%s',[unq{n,1} ' (Median)']);
end

% Here the fold change headings
for n = 1:numel(unq)
    for r = n+1:numel(unq)
        fprintf(fid,'\t%s',['MedFC-' unq{n,1} '-' unq{r,1}]);
    end
end

fprintf(fid,'\n');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [name,class,mzDB,ppmDB] = genAnnoInfo(match,db,dbmz,mz)

numM = numel(match);

for n = 1:numM
    
    % Annotation name
    if isfield(db,'Name')
        tmpName = db.Name{match(n)};
    elseif isfield(db,'Formula')
        tmpName = db.Formula{match(n)};
    else
        tmpName = '?';
    end
    if n == 1
        name = tmpName;
    else
        name = [name ' / ' tmpName];
    end
    
    % Class
    tmpClass = [db.Class1{match(n)} ',' db.Class2{match(n)}];
    if n == 1
        class = tmpClass;
    else
        class = [class ' / ' tmpClass];
    end
    
    % m/z value
    tmpDBMZ = sprintf('%0.4f',dbmz(n));
    if n == 1
        mzDB = tmpDBMZ;
    else
        mzDB = [mzDB ' / ' tmpDBMZ];
    end
    
    % ppm value
    tmpPPM = sprintf('%0.1f',1e6 * (mz - dbmz(n)) / dbmz(n));
    if n == 1
        ppmDB = tmpPPM;
    else
        ppmDB = [ppmDB ' / ' tmpPPM];
    end
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
function [allStat] = determinePQvalues(sp,grp,wb,qThresh)
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
    
    if mod(n,100) == 0
        waitbar(n/numV,wb,[int2str(n) '/' int2str(numV)]);
    end
    
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
        [a,b,c] = anova1(sp(:,n),grp,'off');
        
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
if qThresh > 0.05
    qThresh = 0.05;
end
allStat(:,3) = getBHYqVls(allStat(:,1)',qThresh);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ allAss ] = annotateOP(mz,sp,grp,ass,db,txtP)
%annotateOP - based on annotation2text - convert the annotation results
% into a kind of text file that can be opened in Excel...
%
% INPUTs
% mz    - mz vector
% sp    - spectral matrix of intensities
% grp   - grouping variable for statistics
% ass   - ass structure obtained from annotateMZ.m
% txtP  - somewhere to save the files

% Simple stuff at the beginning. Flag tells is if we do the stats
[wb,fid,flag,txtP] = preamble(grp,txtP);

% pq values
if flag
    [allStat] = determinePQvalues(sp,grp,wb);
end
delete(wb);

% For some of the diagnostics plotting capabilities, we need to save the
% true and measured mz values so we can see about systematic drift (obvs
% assuming that the assignments are in fact correct). So store all the
% things in the cell, then pass to another function for actually doing the
% plotting and stuff.
%allAss = cell(numV,2);

% Determine the group unique-ness
[unq,~,unqIdx] = unique(grp);

% NUmber of adducts
numA = size(ass.match,2);

% Loop through each...
numV = size(sp,2);
for n = 1:numV
    
    if mod(n,100) == 0
        %waitbar(n/numV,wb,[int2str(n) '/' int2str(numV)]);
    end

    % Determine annotation strings
    tmpAnno = ass.match(n,:);
    addString = cell(1,numA);
    for r = 1:numA
        if ~isempty(tmpAnno{r})            
            % Generate the annotation string for this adduct
            [addString{r},annoMZ,annoPPM] = getAddString(lm,tmpAnno,adduct);
        else
            addString{r} = '';
        end
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wb,fid,flag,txtP] = preamble(grp,txtP)

% Wait bar draw
wb = waitbar(0,'Initiating','Name','Exporting Annotations');

% File path - define if empty
if isempty(txtP)
    txtP = [pwd filesep 'LipidAnnotations-' datestr(now,'yymmdd-HHMM') '.txt'];
end
disp(txtP);

% Create the file
fid = fopen(txtP,'w');

% Make function work without stats / spectra
if isempty(grp)
    flag = false;
    fprintf(fid,'\n');
    return
else
    flag = true;
end

% Groups
[unq,~,~] = unique(grp);
if numel(unq) == 2
    statTest = 'MannU pVal';
    statTst2 = 'MannU qVal';
else
    statTest = 'ANOVA pVal';
    statTst2 = 'ANOVA qVal';
end

% Write the column headings...
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s',...
    'm/z','nID','IDs','True m/z','±ppm',...
    'Class','nCarb','nDesat','Iso','Diff',statTest,statTst2);

% Add in the headings for the specific groups...
for n = 1:numel(unq)
    fprintf(fid,'\t%s',[unq{n,1} ' (Mean)']);
end
for n = 1:numel(unq)
    fprintf(fid,'\t%s',[unq{n,1} ' (Median)']);
end

% Here the fold change headings
numFC = sum(1:numel(unq)-1);
for n = 1:numel(unq)
    for r = n+1:numel(unq)
        fprintf(fid,'\t%s',['MedFC-' unq{n,1} '-' unq{r,1}]);
    end
end

fprintf(fid,'\n');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allStat] = determinePQvalues(sp,grp,wb)
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
        [a,b] = anova1(sp(:,n),grp,'off');
        
        % Need to format the output...
        if isnan(b(3))
            allStat(n,1) = NaN;
        else
            allStat(n,1) = b(3);
            
            % Which is the largest group?
            allStat(n,2) = find(strcmp(unq,a) == 1);
            
        end
        
    end
    
end

% Here we do the q-value correction...
allStat(:,3) = getBHYqVls(allStat(:,1)',0.05);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ints] = prepStuff(ass,data,grp,grpIdx,pq,match)
% Here we determine the things that are to be printed

% Now the me(di)an intensities for each of the various groups...
ints = zeros(numel(grp),2);
for n = 1:size(ints,1)
    fx = grpIdx == n;
    ints(n,1) = mean(data(fx));
    ints(n,2) = median(data(fx));
end

nCarb   = [];
nDesat  = [];
trueMZ  = [];
truePPM = [];

% This is the db index
match = ass.match(n,:)

% How many tentative annotations have been made?
numMatches = numel(fx);

% Define these...
if numMatches > 0
    
    
    
else
    % Where there are no annotations, just need blank fields for these
    annoClss = '';
    nCarb = '';
    nDesat = '';
    
end

% Additional things...
annoMZ   = ass.annoMZ(1,n);
annoName = ass.annoNam{n,1};
iso      = [];

% Add the annotiation / assignments into allAss
allAss{n,1} = op.cmz(1,n);
allAss{n,2} = trueMZ;

diffGrp = allStat(n,2);
if diffGrp == 0
    diffGrp = 'NA';
else
    diffGrp = unq{allStat(n,2)};
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [annoClass,annoMZ,annoPPM] = getAddString(lm,match,adduct)
% For each adduct possibility, determine adduct annotation string...

numMatches = numel(match);

for r = 1:numMatches
    
    % Index of the annotated match
    i = match{r};

    % Class...
    tmp = [lm.Class1{i} ',' lm.Class{2} ',' lm.Class{3}];
    if r == 1
        annoClass = tmp;
    else
        annoClass = [annoClass ' / ' tmp];
    end
    
    % m/z
    trueMass = lm.Mass(i);
    trueMZ   = mass2mz(trueMass,adduct);
    tmpPPM   = 1e6 * (measMZ - trueMZ) / trueMZ;
    
    tmpMZ = sprintf('%0.4f',trueMZ);
    tmpPPM = sprintf('%0.1f',tmpPPM);
    if r == 1
        annoMZ = trueMZ;
        annoPPM = truePPM;
    else
        annoMZ  = [annoMZ  ', ' tmpMZ];
        annoPPM = [annoPPM ', ' tmpPPM];
    end
    
    
end

return

for n = 1:2
    % Get the classes of the assigned lipids, and if more than one,
    % append to the end of the existing one...
    if r == 1
        annoClss = [ass.db.c1{dbi(r)} ',' ass.db.c2{dbi(r)}];
    else
        annoClss = [annoClss ' / ' ...
            ass.db.c1{dbi(r)} ',' ass.db.c2{dbi(r)}]; %#ok<*AGROW>
    end

    % Now get the true m/z values of the assignments
    tmpTrue = ass.db.mz(dbi(r));
    tmpPPM  = 1e6 * (op.cmz(1,n) - tmpTrue) / tmpTrue;
    if r == 1
        trueMZ = sprintf('%0.4f',tmpTrue);
        truePPM= sprintf('%0.1f',tmpPPM);
    else
        trueMZ = [trueMZ ', ' sprintf('%0.4f',tmpTrue)];
        truePPM= [truePPM ', ' sprintf('%0.1f',tmpPPM)];
    end



    % Now work out the length/desaturations in the thing
    tmp = ass.db.nam{1,dbi(r)};

    cho = strfind(tmp,'(');
    chc = strfind(tmp,')');
    chl = strfind(tmp,':');

    if ~isempty(cho) && ~isempty(cho) && ~isempty(chl)
        if isempty(nCarb)
            nCarb  = tmp(cho(1)+1:chl(1)-1);
            nDesat = tmp(chl(1)+1:chc(1)-1);
        else
            nCarb  = [nCarb ' // ' tmp(cho(1)+1:chl(1)-1)];
            nDesat = [nDesat ' // ' tmp(chl(1)+1:chc(1)-1)];
        end
    else
        % Declare yourself flummoxed
        if isempty(nCarb)
            nCarb = '(?)';
            nDesat = '(?)';
        else
            nCarb = [nCarb ' // (?)'];
            nDesat = [nDesat ' // (?)'];
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
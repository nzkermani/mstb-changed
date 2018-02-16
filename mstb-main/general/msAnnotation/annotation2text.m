function [ allAss ] = annotation2text(ass,txtP,op)
%annotation2text - convert the annotation results into a kind of text file
%that can be opened in Excel... The ph is a structure for the writing of
%the anova/post hoc analyses out with the files as well. Need to include as
%appropriate...

wb = waitbar(0,'Initiating','Name','Exporting Annotations');

% File path
if isempty(txtP)
    txtP = [pwd filesep 'NimaLipids-' datestr(now,'yymmdd-HHMM') '.txt'];
end
disp(txtP);

% Create the file
fid = fopen(txtP,'w');

% Groups
[unq,~,unqInd] = unique(op.histID);
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

% Convert the sparse op.XPeaks matrix to a full one
op.XPeaks = full(op.XPeaks);

% How many variables are there?
numV = size(op.XPeaks,2);

% For some of the diagnostics plotting capabilities, we need to save the
% true and measured mz values so we can see about systematic drift (obvs
% assuming that the assignments are in fact correct). So store all the
% things in the cell, then pass to another function for actually doing the
% plotting and stuff.
allAss = cell(numV,2);

% In order to do correction of the p values into q values, we have to
% calculate the p-values all together initially...
allStat = zeros(numV,3);
for n = 1:numV
    
    if mod(n,100) == 0
        waitbar(n/numV,wb,[int2str(n) '/' int2str(numV)]);
    end
    
    % Get the statistics for this variable, either MWU or ANOVA (n > 2)
    if numel(unq) == 2
        % Do mannU        
        
        iX = unqInd == 1;
        iY = unqInd == 2;
        
        
        % Get the p value in sv
        [allStat(n,1)] = ranksum(op.XPeaks(iX,n),op.XPeaks(iY,n));
        
        % Determine the group in which it is biggest        
        if median(op.XPeaks(iX,n)) > median(op.XPeaks(iY,n))
            allStat(n,2) = 1;
        else
            allStat(n,2) = 2;
        end
        
    else
        % Do ANOVA
        [a,b] = anovaPH(op.cmz(n),op.XPeaks(:,n),op.histID);
        
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
    

% Loop through each...
for n = 1:numV
    
    if mod(n,100) == 0
        waitbar(n/numV,wb,[int2str(n) '/' int2str(numV)]);
    end
    
    % Everything variable gets printed in the file. We start with the
    % simple things about the potential annotations (if there is/are
    % one/any), then the results of the anovaPH and finally the mean/median
    % values for each annotation...
    
    % EVERY SINGLE VARIABLE MUST BE PRINTED!
    
    % Determine the group in which the variable is largest...
%     if numel(unq) == 2
%         % Do mannU        
%         
%         iX = unqInd == 1;
%         iY = unqInd == 2;
% %         
% %         
% %         % Get the p value in sv
% %         [sv] = ranksum(op.XPeaks(iX,n),op.XPeaks(iY,n));
% %         sv = [0 0 sv];
%         
%         % Determine the group in which it is biggest        
%         if median(op.XPeaks(iX,n)) > median(op.XPeaks(iY,n))
%             sn = unq(1);
%         else
%             sn = unq(2);
%         end
%         
%         
%         
%     else
%         % Do ANOVA
%         [sn,sv] = anovaPH(op.cmz(n),op.XPeaks(:,n),op.histID);
%     end
        
    % Don't want to print pvalues equal to 0
    if allStat(n,1) == 0
        pval = [];
    elseif isnan(allStat(n,1))
        pval = [];
    else
        pval = allStat(n,1);
    end
    
    % Now the me(di)an intensities for each of the various groups...
    ints = zeros(numel(unq),2);    
    for r = 1:size(ints,1)
        fx = strcmp(op.histID,unq{r,1});
        fx = find(fx == 1);
        ints(r,1) = mean(op.XPeaks(fx,n));
        ints(r,2) = median(op.XPeaks(fx,n));
    end
    
    nCarb   = [];
    nDesat  = [];
    trueMZ  = [];
    truePPM = [];
    
    % This is the db index
    [fx,~] = find(ass.list(:,1) == n);
    
    % How many tentative annotations have been made?
    numMatches = numel(fx);
                
    % Define these...
    if numMatches > 0
        
        dbi = ass.list(fx,2);
        for r = 1:numMatches
            
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
    
    % Finally print the stuff
    fprintf(fid,'%0.4f\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%0.6E\t%0.6E',...
        annoMZ,numMatches,annoName,...
        trueMZ,truePPM,...
        annoClss,deblank(nCarb),deblank(nDesat),...
        iso,diffGrp,pval,allStat(n,3));

    % Add in the mean intensities...
    for s = 1:numel(unq)
        fprintf(fid,'\t%0.4f',ints(s,1));
    end
    for s = 1:numel(unq)
        fprintf(fid,'\t%0.4f',ints(s,2));
    end
    
    % Here print out the ME++DIAN fold changes
    for n = 1:numel(unq)
        for r = n+1:numel(unq)
            tmpFC = log2(ints(n,2) / ints(r,2));
            fprintf(fid,'\t%0.4f',tmpFC);            
            %fprintf(fid,'\t%s',['MedFC-' unq{n,1} '-' unq{r,1}]);
        end
    end

    % New line
    fprintf(fid,'\n');

    
    
    
end

% Close the text file
fclose(fid);

delete(wb);

end


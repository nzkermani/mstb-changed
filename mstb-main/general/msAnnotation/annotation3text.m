function [ allAss ] = annotation3text(ass,txtP,mz)
%annotation3text - toned down version of its predecessor

if isstruct(mz)
    mz = mz.cmz;
end

wb = waitbar(0,'Initiating','Name','Exporting Annotations');

% Say to where you are exporting the data
disp(txtP);

% Create the file
fid = fopen([txtP '.txt'],'w');


% Write the column headings...
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s',...
    'm/z','nID','IDs','True m/z','±ppm','Class','nCarb','nDesat');
fprintf(fid,'\n');

% How many variables are there?
numV = numel(mz);

% For some of the diagnostics plotting capabilities, we need to save the
% true and measured mz values so we can see about systematic drift (obvs
% assuming that the assignments are in fact correct). So store all the
% things in the cell, then pass to another function for actually doing the
% plotting and stuff.
allAss = cell(numV,2);

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
    clear annoMZ numMatches annoName annoClss
    
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
            tmpPPM  = 1e6 * (mz(n) - tmpTrue) / tmpTrue;
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
    annoMZ   = ass.annoMZ(n);
    annoName = ass.annoNam{n,1};
    iso      = [];

    % Add the annotiation / assignments into allAss
    allAss{n,1} = mz(n);
    allAss{n,2} = trueMZ;
    
    % Finally print the stuff
    fprintf(fid,'%0.4f\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d',...
        annoMZ,numMatches,annoName,...
        trueMZ,truePPM,...
        annoClss,deblank(nCarb),deblank(nDesat));

    % New line
    fprintf(fid,'\n');

    
    
    
end

% Close the text file
fclose(fid);

delete(wb);

end


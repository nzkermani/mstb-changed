function statsUnivariateCalculate(~,~,fig,window)
%statsLDACalculate - needs to be separate from the drawing function, as LDA
%changes with different groupings...


% Guidata
sts = guidata(fig.fig);

% Determine the type of analysis
val = window.analysis.Value;
str = window.analysis.String;
analysis = str{val};

% Determine the groups/groups that are to be differentiated
val = window.groups.Value;
str = window.groups.String;
grp = statsObservationLabels(sts.proc.meta,str,val);
numG = numel(unique(grp));

% There must be at least two groups
if numG < 2
    disp('Cannot calculate with less than 2 groups');
    return
elseif strcmp(analysis,'Output') && numG ~= 2
    disp('Output only possible with 2 groups');
    return   
end
        
% Need to extract the pixels from the data, run ANOVA on them and
% find a way to plot the data...
switch analysis
    case 'Text Output'
        
        % Warn that in the case of log transformed data, the results may
        % look a little off, as the lowest value is no longer zero and it
        % may not give the best results
        if sts.proc.opts.doLog
            msg = ['The data is log transformed. '...
                'You may wish to un-transform the data and start '...
                'again.  If you wish to proceed, click accordingly'];
            
            choice = questdlg(msg,'Question?',...
                'Cancel','Proceed Anyway',...
                'Cancel');
            
            if strcmp(choice,'Cancel')
                return
            end
        end
        
        annotate(sts.proc.var.mz,sts.proc.sp,grp);
        
        return
        
    otherwise
        [pq,fc] = statsUnivariateMethod(sts.proc.var.mz,sts.proc.sp,...
            grp,analysis);
end

% Format results for output
sts.res.uv.method = analysis;
sts.res.uv.grp = grp;
sts.res.uv.pq = pq;
sts.res.uv.fc = fc';

% Save the results
guidata(fig.fig,sts);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function export2text(window)

% This part is taken from the matlab save function and just allows us to
% use the -v7.3 flag as default for large file sizes...
filters = {'*.txt','txt files (*.txt)'};

% Determine default file names etc...
defP = [pwd filesep];
defN = ['Annotations-' datestr(now,'yymmdd-HHMMSS')];

% Just ask the user for a file name and location...
[fn,pn,~] = uiputfile(filters,...
    getString(message('MATLAB:uistring:filedialogs:SaveWorkspaceVariables')),...
    [defP defN '.txt']);

% If cancelled then nothing happens...
if isequal(fn,0)
    return
end

% Now continue with the function...
pqThresh = str2double(get(window.pqThresh,'String'));

% Get the database - this needs to be specified as we don't know what
[db,ass] = annotateMZ(mz,...
    'Polarity','negative',...
    'Adduct',{'M-H','M+Cl','M-H2O-H'},...
    'Tolerance',8,...
    'Database','dipa');

% Output to a text file
annotateOP2(mz,sp,...
    histID,...
    ass,...
    db,...
    txtP,...
    qThresh);




    
    % Create the full length file name
    fn = strrep(fullfile(pn,fn), '''', '''''');
    disp(fn)
    
    % Save the file name - need to test that 'normal' files are ok with
    % this too.
    %save(fn,'dpn','-v7.3');
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

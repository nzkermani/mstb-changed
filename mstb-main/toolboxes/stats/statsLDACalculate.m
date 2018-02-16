function statsLDACalculate(~,~,fig,window)
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

% Simple check to ensure that we don't waste our time with too many groups,
% or perhaps where there are too few observations per group
numG = numel(unique(grp));
if numG > 30 || numG == 1 || numG == size(grp,1)
    disp('Too many / too few groups for supervised analysis');
    wd = warndlg('Too many / too few groups for analysis');
    pause(1)
    delete(wd);
    return
end

% Create a waitbar
wb = waitbar(0.25,'LDA etc...');

% Now run the analyses...
switch analysis
    
    case 'MMC-LOO'
        
        % Run MMC with LOOCV
        [ss,cm,preds,grps,ll] = looMMC(sts.proc.sp,grp);
        
        % Save the various elements...
        sts.res.mmc.ss = ss;
        sts.res.mmc.ll = ll;
        sts.res.mmc.cm.cm = cm;
        sts.res.mmc.cm.names = grps;
        sts.res.mmc.grp = grp;
        
    case 'MMC-k10'
        
        % Run MMC with 10-fold CV
        [ss,cm,preds,grps,ll] = kfoldMMC(sts.proc.sp,grp,10);
        
        % Save the various elements...
        sts.res.mmc.ss = ss;
        sts.res.mmc.ll = ll;
        sts.res.mmc.cm.cm = cm;
        sts.res.mmc.cm.names = grps;
        sts.res.mmc.grp = grp;
        
    case 'MMC-LPO'
        
        % Determine the group to be used for leave <patient> out
        val = window.cvgroup.Value;
        str = window.cvgroup.String;
        cvg = sts.proc.meta.(str{val});
        
        % Now here we run the analysis
        [ss,cm,preds,grps,ll,roc] = lpoMMC(sts.proc.sp,grp,cvg);
        
        assignin('base','roc',roc);
        if ~isempty(roc)
            set(window.roc,'Visible','on');
        end
        
        % Save the various elements...
        sts.res.mmc.ss = ss;
        sts.res.mmc.ll = ll;
        sts.res.mmc.cm.cm = cm;
        sts.res.mmc.cm.names = grps;
        sts.res.mmc.grp = grp;
        sts.res.mmc.roc = roc;
                
    case 'MMC-Ext'
        
        % This is for externally validated MMC, so run it with no internal
        % cross validation. Instead, we expect that this is performed
        % elsewhere with different samples. Somehow it will be possible
        % here...
        numDV = max([numG-1 2]);
        [~,~,unqInd] = unique(grp);
        [~,ss,~,ll] = recursiveMmcLda(sts.proc.sp,unqInd,numDV);
        
        % Save the various elements
        sts.res.mmc.ss = ss;
        sts.res.mmc.ll = ll;
        sts.res.mmc.cm = [];
        sts.res.mmc.grp = grp;

        
end

% Set the maxComp value
maxComp = size(ss,2);
try
    set(window.comp1,'Value',1,'String',int2str([1:maxComp]'));
    set(window.comp2,'Value',2,'String',int2str([1:maxComp]'));
catch
end

% Update GUI
guidata(fig.fig,sts);

% Here is where we trigger the drawing function
statsLDADraw([],[],fig,window);

% Delete the waitbar, although it appears to go before the plot has been
% properly updated
delete(wb);

end


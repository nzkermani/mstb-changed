function desiStatisticsDraw(~,event,fig,man)
% desiStatisticsDraw in the window...

% Guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Check that there is actucally something to statisticify
if ~isfield(dpn,'anno')
    disp('No annotations');
    return
end

% Determine the option selected
val = get(man.list,'Value');
str = get(man.list,'String');
method = str{val};

% If we are doing the PCA/MMC of all unambiguously classified pixels, we
% need to consider this.
allPixels = ~isempty(strfind(method,'Unambig.'));

% Get the necessary pixels/identities
[sp,histID,~,~] = extractAnnotatedPixels(dpn,man,allPixels,'d1');
if isempty(sp)
    disp('NO GROUPS SELECTED FOR STATISTICS');
    return
end

% What about the second set of data if doing dual mode?
if strcmp(dpn.mode,'dual')
    [sp2,~,~,mask] = extractAnnotatedPixels(dpn,man,allPixels,'d2');
    [~,~,isInt1,isInt2,pixID] = dpnAnnotationExtract(dpn);
    isInt1 = isInt1(mask > 0,:);
    isInt2 = isInt2(mask > 0,:);
    pixID  = pixID(mask  > 0,:);
    
end

% Now is the time to plot stuff...
if strcmp(get(fig.tb.layout,'State'),'off')
    set(fig.tb.layout,'State','on');
    if strcmp(dpn.mode,'single')
        desiChangeLayout(fig.tb.layout,event,fig);
    else
        dpnChangeLayout(fig.tb.layout,event,fig);
    end
end

% Remove any button down function clicks from the axes
set(fig.ax.mv,'ButtonDownFcn','',...
    'XTick',[],...
    'XTickLabel',[],...
    'YTick',[],...
    'YTickLabel',[]);
xlabel(fig.ax.mv,'');
set(fig.ax.sp,'ButtonDownFcn','',...
    'XTick',[],...
    'XTickLabel',[],...
    'YTick',[],...
    'YTickLabel',[]);
xlabel(fig.ax.sp,'');

% Here we need to do the stuff related to fusion. Everything in one go to
% make life a little easier.
if strcmp(dpn.mode,'dual')   
    
    % Determine the desired method
    fv = get(man.fuse,'Value');
    ft = get(man.fuse,'String');
    fuseMethod = ft{fv};
    
    % Run the function for LL fusion - note the HL comes later on
    if strcmpi(fuseMethod,'ll (concat)');
        [mzFuse,spFuse] = prepFusion(dpn.d1.mz,sp,...
            dpn.d2.mz,sp2,man,fuseMethod);   
    end
    
end


% Enforce non-transformation of the data for univariate analyses
if strcmpi(method,'anova') || ...
        strcmpi(method,'kruskal-wallis') || ...
        strcmpi(method,'annotate') || ...
        strcmpi(method,'spectra----')
    set(man.log,'Value',1);
    set(man.scale','Value',1);    
else
    set(man.volcano,'Value',0);
end
doUVvolcano = get(man.volcano,'Value');

% Here we implement the norm/tran/scal procedures
[mz,sp] = xxxNormTranScal(man,dpn.d1.mz,sp);
%[mz,sp] = xxxTranNormScal(man,dpn.d1.mz,sp);
if strcmp(dpn.mode,'dual')
    [mz2,sp2] = xxxNormTranScal(man,dpn.d2.mz,sp2);
    %[mz2,sp2] = xxxTranNormScal(man,dpn.d2.mz,sp2);
end

% What about high level fusion?
retroLoad = false;
retroScores = false;
if strcmp(dpn.mode,'dual')
    if strcmpi(fuseMethod(1:2),'hl')        
        
        [mzFuse,spFuse,hlFuse] = prepFusion(dpn.d1.mz,sp,...
            dpn.d2.mz,sp2,man,fuseMethod);
        
        % What is the choice regarding loadings plots?
        retroLoad = true(get(man.retroLoadings,'Value'));
        
        retroScores = true(get(man.retroScores,'Value'));
        
    end
end

% Determine stuff for the univariate options - set these to 0 if not doing
% univariate analysis


% Main plot
switch method
    
    case 'Export'
        
        % Dump the annotated pixels to the workspace
        export.mz = mz;
        export.sp = sp;
        export.meta.histID = histID;
        
        assignin('base','desiExport',export);
        return
        
    case {'PCA','MMC','PCA-Unambig.','MMC-Unambig.','MMC-LOO','MMC-k10'}
        
        % Either single or dual mode analysis?
        if strcmp(dpn.mode,'dual')
            
            % Run the statistics part
            cs = [1 2];
            [llF,ssF,~,ftS,~] = multivariateMethod(man,spFuse,histID,method,cs);
            
            % Save the scores
            dpn.fuse.mva.(method).sc = ssF;
            
            
            
            % What about retroscores? The aim is perhaps to calculate the
            % scores of the separate datasets using the loadings...
            if retroScores
                
                % Recalculation of scores
                tmpSS1 = bsxfun(@minus,sp, mean(sp,1))  * hlFuse(1).ll * llF(1:size(hlFuse(1).ss,2),:);
                tmpSS2 = bsxfun(@minus,sp2,mean(sp2,1)) * hlFuse(2).ll * llF(size(hlFuse(1).ss,2)+1:end,:);

                % Dual scatter plot
                doScatterPlotDual(fig.ax.mv,dpn.anno,...
                    tmpSS1,tmpSS2,ftS,histID);
                
            else
                % Scores - just the fused ones now
                doScatterPlot(fig.ax.mv,dpn.anno,ssF,ftS,histID,cs);
            end
            
            % Loadings - could do retro projection of loadings or just show
            % the loadings of the high level PCA scores - it all depends on
            % the method selected...
            if retroLoad
                
                % Recalculation of loadings
                hlFuse(1).ll = hlFuse(1).ll * llF(1:size(hlFuse(1).ss,2),:);
                hlFuse(2).ll = hlFuse(2).ll * llF(size(hlFuse(1).ss,2)+1:end,:);
                
                % Combine into single vector...
                llF = vertcat(hlFuse.ll);
                mzFuse = [mz mz2];
                
                doLoadingsPlotDual(fig.ax.sp,mzFuse,llF,fuseMethod,histID);
                
                % Save the high-level loadings
                dpn.fuse.mva.(method).ld.mz1 = mz;
                dpn.fuse.mva.(method).ld.mz2 = mz2;
                dpn.fuse.mva.(method).ld.d1 = hlFuse(1).ll;
                dpn.fuse.mva.(method).ld.d2 = hlFuse(2).ll;
                

            else
                % Just show the loadings of the second set of scores
                doLoadingsPlotDual(fig.ax.sp,mzFuse,llF,fuseMethod,histID);            
                
                % Save the low-level loadings
                cu = numel(mz);
                dpn.fuse.mva.(method).ld.mz1 = mz;
                dpn.fuse.mva.(method).ld.mz2 = mz2;
                dpn.fuse.mva.(method).ld.d1 = llF(1:cu,:);
                dpn.fuse.mva.(method).ld.d2 = llF(cu+1:end,:);

            end
            
        else
            
            % Run the statistics part
            cs = [1 2];
            [ll,ss,ee,titScore,titLoad] = multivariateMethod(man,sp,histID,method,cs);

            % Scores
            doScatterPlot(fig.ax.mv,dpn.anno,ss,titScore,histID,cs);

            % Loadings
            doLoadingsPlot(fig.ax.sp,mz,ll,titLoad,histID);
                
            % Make me a confusion matrix if we have done MMC
            if isstruct(ee)
                % Go!
                
                % Need to determine the colours
                [unq,~,ind] = unique(histID);
                numG = numel(unq);
                cols = zeros(numG,3);
                for n = 1:numG

                    % Indices of points to plot
                    fx = ind == n;

                    % Colour finding/matching
                    cx = strcmp(dpn.anno(:,5),unq{n});
                    cx = find(cx == 1,1,'first');
                    cols(n,:) = dpn.anno{cx,3};
                end

                
                % Mark incorrect predictions on the plot
                fx = ee.pred(:,1) ~= ee.pred(:,2);
                scatter(fig.ax.mv,ss(fx,1),ss(fx,2),250,'k','o');
                
                % What about confusion matrices?
                doConfMat = get(man.confmat,'Value');
                if doConfMat
                    makeConfMat([],ee.confMat,cols,[],ee.groups);
                    makeConfMat([],ee.confMat,cols,[],[]);
                end
            end
        end
        

    case 'Spectra'
        doSpectralPlot(fig.ax.sp,mz,sp,'Spectra',histID,dpn.anno);
        doHeatMap(fig.ax.mv,mz,sp,'Heat Map',histID,dpn.anno);
        
    case {'ANOVA','Kruskal-Wallis'}
        
        % Need to extract the pixels from the data, run ANOVA on them and
        % find a way to plot the data...
        [pq,fc] = univariateMethod(man,mz,sp,histID,method);
        
        % Find out if we wanted a volcano or mz plot
        if doUVvolcano
            doVolcanoPlot(fig.ax.mv,...
                fig.ax.sp,...
                mz,fc',pq(:,2),...
                [method 'Volcano Plot of log2 fold changes against -log10 q-values'],...
                histID,sp,dpn.anno);
        else        
            doStemPlot(fig.ax.mv,...
                fig.ax.sp,...
                mz,fc',pq(:,2),...
                [method ' Stem Plot of log2 fold changes against -log10 q-values'],...
                histID,sp,dpn.anno);
        end
        
    case 'Annotate'
        % Function to export the stuff with annotations, p/q values and a
        % bunch of other stuff...
        disp('Performing annotations');
        
        % Determine the normalisation method
        val2 = get(man.norm,'Value');
        lab = get(man.norm,'String');
        normMethod = lab{val2};
        
        % Annotate window instead
        annotate(mz,sp,histID,...
            'Unicorn',rand(1) < 0.5,...
            'preNorm',normMethod,...
            'Folder',dpn.file.dir,...
            'File',datestr(now,'yymmdd-HHMMSS'));
        
    case {'PCA-IntPos','PCA-IntNeg'}
        
        % Run the statistics part
        if strcmp(method,'PCA-IntPos')
            useMZ = mz;
            useData = sp;
            c2one = classMany2One([histID isInt1]);
        else
            useMZ = mz2;
            useData = sp2;
            c2one = classMany2One([histID isInt2]);
        end
        
        cs = [1 2];
        [ll,ss,ee,titScore,titLoad] = multivariateMethod(man,useData,c2one,'PCA',cs);

        % Scores
        doScatterPlot(fig.ax.mv,dpn.anno,ss,titScore,c2one,cs);

        % Loadings
        doLoadingsPlot(fig.ax.sp,useMZ,ll,titLoad,histID);
        
    case {'MMC-IntPos','MMC-IntNeg'}
                
        % Run the statistics part
        delim = ', ';
        if strcmp(method,'MMC-IntPos')
            useMZ = mz;
            useData = sp;
            c2one = classMany2One([histID isInt1],delim);
        else
            useMZ = mz2;
            useData = sp2;
            c2one = classMany2One([histID isInt2],delim);
        end

        % Run the standard multivariate method
        [ll,ss,ee,titScore,titLoad] = multivariateMethod(man,useData,c2one,'MMC');
        
        % Scores
        cs = [1 2];
        [grpInfo] = doScatterPlot(fig.ax.mv,dpn.anno,ss,titScore,c2one,cs,delim);

        % Loadings
        doLoadingsPlot(fig.ax.sp,useMZ,ll,titLoad,histID);

        % Mark incorrect predictions on the plot
        %fx = ee.pred(:,1) ~= ee.pred(:,2);
        %scatter(fig.ax.mv,ss(fx,1),ss(fx,2),250,'k','o');

        % What about confusion matrices?
        doConfMat = get(man.confmat,'Value');
        
        if doConfMat
            makeConfMat([],ee.confMat,grpInfo,[],ee.groups);
            makeConfMat([],ee.confMat,grpInfo,[],[]);
        end
        
    case {'Annotation-Stats'}
        % Need to run some function to determine information about the
        % annotation statistics... Have to base this on annotations, those
        % pixels, and ...
        desiAnnotationStats(dpn);
        
    otherwise
        error('There is no otherwise');
end

% Save the guidata in case there were various updates
guidata(fig.fig,dpn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mzF,spF,hl] = prepFusion(mz1,sp1,mz2,sp2,man,method)
% Get all the data required for data fusion

switch lower(method(1:2))
    case 'll'
        
        % Concatenate the two matrices together (only ever two)
        mzF = [mz1 mz2];
        spF = [sp1 sp2];
        
        % Perform the selected norm/tran/scal methods
        [mzF,spF] = xxxNormTranScal(man,mzF,spF);       
        
        hl = [];
    case 'hl'
                
        % Determine if we need to select on the number of components or the
        % percentage of variance expressed
        if strcmpi(method(4:end),'(pca #)')
            maxComp = str2double(get(man.fuseComp,'String'));
            
        elseif strcmpi(method(4:end),'(pca %)')
            maxComp = size(sp1,1);
            
        end
        
        
        % Both of sp2 and sp2 have already been norm/tran/scal in the same
        % method. Thus all we need to do is run PCA, then combine xPCs and
        % scale according to variance. Also calculate a PC1-derived
        % variance scaling factor, and divide the scores and loadings by it
        [hl(1).ll,hl(1).ss,hl(1).ee] = pca(sp1,'NumComponents',maxComp);
        hl(1).ee = 100 * hl(1).ee / sum(hl(1).ee);
        hl(1).vr = std(hl(1).ss(:,1),[],1);
        hl(1).ss = hl(1).ss / hl(1).vr;
        hl(1).ll = hl(1).ll / hl(1).vr;
        
        [hl(2).ll,hl(2).ss,hl(2).ee] = pca(sp2,'NumComponents',maxComp);
        hl(2).ee = 100 * hl(2).ee / sum(hl(2).ee);
        hl(2).vr = std(hl(2).ss(:,1),[],1);
        hl(2).ss = hl(2).ss / hl(2).vr;
        hl(2).ll = hl(2).ll / hl(2).vr;
        
        % Here is where we need to decide on how many components to include
        % if doing it based on the % of variance
        if maxComp == size(sp1,1)
            maxPerc = str2double(get(man.fuseComp,'String'));
            
            numComp(1) = min([numel(hl(1).ee) find(cumsum(hl(1).ee) > maxPerc,1,'first')]);
            numComp(2) = min([numel(hl(2).ee) find(cumsum(hl(2).ee) > maxPerc,1,'first')]);
            
            if numComp(1) == 1
                numComp(1) = 2;
            end
            if numComp(2) == 1
                numComp(2) = 2;
            end
            
            % Trim out the crap
            for n = 1:2
                hl(n).ll = hl(n).ll(:,1:numComp(n));
                hl(n).ss = hl(n).ss(:,1:numComp(n));
                hl(n).ee = hl(n).ee(1:numComp(n));
            end
            
        else
            numComp = [maxComp maxComp];            
        end
        
        % Now we have two datasets with appropriately* scaled PC scores.
        % These just need to be combined...
        mzF = [1:numComp(1) 1:numComp(2)];
        spF = horzcat(hl.ss);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp,histID,pixID,mask2] = extractAnnotatedPixels(dpn,man,allPixels,flag)
% Get the pixels that were annotated and requested to be included in this
% analysis

% Determine which annotated groups are to be included according to the
% table
tabdat = get(man.tab,'Data');
check = cellfun(@eq,tabdat(:,1),repmat({1},[size(tabdat,1) 1]));
if sum(check) == 0
    sp = [];
    histID = [];
    pixID = [];
    return
end

% Extract the annotated pixels from the data...
[mask2,histID,pixID] = desiAnnotationExtract(dpn);

% Decide if we need to plot all unambiguous pixels rather than just the
% annotated ones.
if allPixels
    try
        mask2 = dpn.d1.mva.MMC.mmcMask2;
        histID = dpn.d1.mva.MMC.mmcClass;
    catch err
        allPixels = false;
    end
end

% Set pixels in mask2 to zero if they weren't selected in the table
numG = size(tabdat,1);
for n = 1:numG
            
    % Find the number from within the string...
    f1 = strfind(tabdat{n,2},'<TD>') + 4;
    f2 = strfind(tabdat{n,2},'</TD>') - 1;
    fx = str2double(tabdat{n,2}(f1:f2));

    % Find unselected groups
    if check(n,1) == 0        
        
        % Find all values of this
        fy = mask2 == fx;
        
        % Set all values of this annotation to zero
        mask2(fy,:) = 0;
        
    else
        
        % Ensure that all have the same annotation
        fy = mask2 == fx;
        histID(fy,:) = tabdat(n,3);
        
        
    end
    
    
end

% Extract the MS data...
switch flag
    case 'd1'        
        sz = size(dpn.d1.sp);
        sp = reshape(dpn.d1.sp,[sz(1)*sz(2) sz(3)]);
        
    case 'd2'
        sz = size(dpn.d2.sp);
        sp = reshape(dpn.d2.sp,[sz(1)*sz(2) sz(3)]);        
end
sp = sp(mask2 > 0,:);

% Generate a cell vector of histological ID / class, which is based on the
% text annotations in dpn.anno
histID = histID(mask2 > 0,:);
pixID = pixID(mask2 > 0,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ll,ss,ee,titText,titLoad] = multivariateMethod(man,sp,histID,method,cs)

wb = waitbar(1/4,'Performing multivariate analysis');

% What can we do?
switch method
    
    case {'PCA','PCA-Unambig.'}
        
        [ll,ss,ee] = pca(sp,'Economy',true);
        ee = 100 * ee / sum(ee);
        titText = ['PCA: PC' int2str(cs(1)) ' (' sprintf('%0.1f',ee(cs(1))) ...
            '%) v PC' int2str(cs(2)) ' (' sprintf('%0.1f',ee(cs(2))) '%)'];
        titLoad = ['Loadings for PC' int2str(cs(1)) ' (top) and PC' int2str(cs(2)) ' (bottom)'];  

    case {'MMC','MMC-LOO'}
        
        % This performs leave one out CV MMC
        [ss,cm,preds,grps,ll] = looMMC(sp,histID);
        
        % Save the confused elements in a structure
        ee.groups = grps;
        ee.pred = preds;
        ee.confMat = cm;
        
        % This is the default text
        titText = 'MMC: Cross validated scores for LV1 v LV2';
        titLoad = 'MMC Weights for LV1 (top) and LV2 (bottom)';
    
    case {'MMC-k10'}
        
        % This performs leave one out CV MMC
        [ss,cm,preds,grps,ll] = kfoldMMC(sp,histID,10);
        
        % Save the confused elements in a structure
        ee.groups = grps;
        ee.pred = preds;
        ee.confMat = cm;
        
        % This is the default text
        titText = 'MMC: Cross validated scores for LV1 v LV2';
        titLoad = 'MMC Weights for LV1 (top) and LV2 (bottom)';
                
        
    case {'MMC-NoCV','MMC-Unambig.'}
                
        % This is the simple MMC, with no CV.
        [grps,~,grpIdx] = unique(histID);
        [~,ss,ee,ll] = recursiveMmcLda(sp,grpIdx,max([2 numel(grps)-1]));

        titText = 'MMC: Scores for LV1 v LV2';
        titLoad = 'MMC Weights for LV1 (top) and LV2 (bottom)';

    otherwise
        % There is no otherwise!
        ll = [];
        ss = [];
        ee = [];
        return
end

delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pq,hi] = univariateMethod(man,mz,sp,histID,method)
% Calculate univariate statistics on the data

% Somewhere to store the results
numV = size(sp,2);
pq = zeros(numV,2);

wb = waitbar(0,'Univariate calculations!');

warning off all

% Run the various methods
switch method    
    case {'ANOVA'}
        for n = 1:numV
            pq(n,1) = anova1(sp(:,n),histID,'off');
            waitbar(n/numV,wb);
        end
        
    case {'Kruskal-Wallis'}
        for n = 1:numV
            pq(n,1) = kruskalwallis(sp(:,n),histID,'off');
            waitbar(n/numV,wb);
        end
        
    otherwise
        error('there is no otherwise');
end

warning on all

% Perform FDR on the p values
[pq(:,2),~] = getBHYqVls(pq(:,1)',0.05);

% Whilst there are possibly multiple groups in this analysis, let's try to
% calculate a fold change of some description
[unq,~,ind] = unique(histID);
numG = numel(unq);
fc = zeros(numG,numV);
for n = 1:numG
    
    % Mean spectrum of reference group
    fx = ind == n;
    ms = nanmean(sp(fx,:),1);
    
    % Store the FCs here
    tmpFC = NaN(numG,numV);
    
    for r = 1:numG        
        if n ~= r
            
            % Mean spectrum of this group
            fy = ind == r;
            rf = nanmean(sp(fy,:),1);
            
            % Calcualte the log2 fold changes
            tmpFC(r,:) = log2(ms ./ rf);
        end
    end
    
    % Remove NaN / Inf
    tmpFC(isnan(tmpFC)) = 0;
    tmpFC(isinf(tmpFC)) = 0;
    
    % Find largest/smallest
    lo = nanmin(tmpFC,[],1);
    hi = nanmax(tmpFC,[],1);    
    fx = abs(lo) > hi;    
    hi(fx) = lo(fx);
    
    % These are the biggest FCs cf this group
    fc(n,:) = hi;
    
    
end
       
lo = nanmin(tmpFC,[],1);
hi = nanmax(tmpFC,[],1);
fx = abs(lo) > hi;
hi(fx) = lo(fx);
           
% Delete the waitbar
delete(wb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grpInfo] = doScatterPlot(parent,anno,scores,titText,histID,cs,delim)
% SCatter plot the MVA results

% This is a simple fudge
if nargin == 6
    delim = '@';
end

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);

% Hold the axes
axes(parent);
hold on;

% Determine unique groups in the histIDs
[unq,~,ind] = unique(histID);
numG = numel(unq);
grpInfo = cell(numG,2);
for n = 1:numG
    
    % Indices of points to plot
    fx = ind == n;
    
    % Colour finding/matching
    findAt = strfind(unq{n},delim);
    if ~isempty(findAt)
        cx = strcmp(anno(:,5),unq{n}(1:findAt(1)-1));
        if strcmp(unq{n}(findAt(1)+length(delim):end),'Measured')
            symb = 'o';
        else
            symb = 'd';
        end
    else
        cx = strcmp(anno(:,5),unq{n});
        symb = 'o';
    end
    cx = find(cx == 1,1,'first');
    col = anno{cx,3};
    
    grpInfo{n,1} = col;
    grpInfo{n,2} = symb;
    
    % Scatter away
    scatter(scores(fx,cs(1)),scores(fx,cs(2)),120,col,symb,'filled',...
        'MarkerEdgeColor','k');
    
    % Determine the ellipse
    [ell,xy,ab] = error_ellipse(scores(fx,cs),95);
    plot(ell(:,1),ell(:,2),'Color',col,'LineWidth',2);
end

% Add a title...
title(titText);
    
box on;

set(gca,'XTickLabel',[],'YTickLabel',[],...
    'XTick',0,...
    'YTick',0,...
    'LineWidth',5,...
    'TickLength',[0 0]);
grid off;

% Draw lines along the origin
xlim('auto');
ylim('auto');
xl = xlim(gca);
yl = ylim(gca);
line([xl(1) xl(2)],[0 0],'LineStyle',':','Color','k');
line([0 0],[yl(1) yl(2)],'LineStyle',':','Color','k');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doScatterPlotDual(parent,anno,scores,scores2,titText,histID)
% SCatter plot the MVA results

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);

% Hold the axes
axes(parent);
hold on;

% Determine unique groups in the histIDs
[unq,~,ind] = unique(histID);
numG = numel(unq);
for n = 1:numG
    
    % Indices of points to plot
    fx = ind == n;
    
    % Colour finding/matching
    cx = strcmp(anno(:,5),unq{n});
    cx = find(cx == 1,1,'first');
    col = anno{cx,3};
    
    % Scatter away
    scatter(scores(fx,1),scores(fx,2),120,col,'o','filled',...
        'MarkerEdgeColor','k');
    scatter(scores(fx,1),scores(fx,2),120,'k','+');
    
    % Determine the ellipse
    [ell,~,~] = error_ellipse(scores(fx,1:2),95);
    plot(ell(:,1),ell(:,2),'Color',col,'LineWidth',2);
    
    % And for the second lot...
    scatter(scores2(fx,1),scores2(fx,2),120,col,'d','filled',...
        'MarkerEdgeColor','k');
    
    % Determine the ellipse
    [ell,~,~] = error_ellipse(scores2(fx,1:2),95);
    plot(ell(:,1),ell(:,2),'Color',col,'LineWidth',2);
    
    % Determine the super ellipse?
    [sull,~,~] = error_ellipse([scores(fx,1:2); scores2(fx,1:2)],95);
    plot(sull(:,1),sull(:,2),'Color',col,'LineWidth',2);
    
end

% Add a title...
title(titText);
    
box on;

set(gca,'XTickLabel',[],'YTickLabel',[],...
    'XTick',0,...
    'YTick',0,...
    'LineWidth',5,...
    'TickLength',[0 0]);
grid off;

% Draw lines along the origin
xlim('auto');
ylim('auto');
xl = xlim(gca);
yl = ylim(gca);
line([xl(1) xl(2)],[0 0],'LineStyle',':','Color','k');
line([0 0],[yl(1) yl(2)],'LineStyle',':','Color','k');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doVolcanoPlot(parent,parent2,mz,fc,pq,titText,histID,sp,anno)
% SCatter plot the MVA results

% Significant thresholds
pqThresh = -log10(0.01);
fcThresh = 2;

% Log the pq values and use these throughout.  Need to remove the Inf
% values, or set them to be above the others by a slight way
lowPQ = min(pq(pq > 0));
idxLo = pq == 0;
pq(isnan(pq)) = 1;
pq = -log10(pq);
pq(idxLo) = -log10(lowPQ) * 1.2;

% Do below but for the plot above
f0 = get(parent2,'Children');
title(parent2,'');
delete(f0);

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);

% Hold the axes
axes(parent);
hold on;

% Highlight only significant ones
fx = abs(fc) > fcThresh & pq > pqThresh;

if sum(fx) == 0
    fx = false(size(pq));
end

% Scatter away - bad variables
scatter(fc(~fx),pq(~fx),60,[0.2627 0.5765 0.7647],'o','filled',...
    'MarkerEdgeColor',[0.7 0.7 0.7]);
    
% These are the interesting ones
scatter(fc(fx),pq(fx),120,[0.8392 0.3765 0.3020],'o','filled',...
    'MarkerEdgeColor',[0.7 0.7 0.7],...
    'HitTest','on',...
    'ButtonDownFcn',{@volcanoClick,...
    [parent parent2],fc(fx),pq(fx),sp(:,fx),histID,mz(fx),anno(:,[3 5])});

% Add text labels so that we can make sure that the boxplots are shown
% correctly
%text(fc(fx),pq(fx),num2str(mz(fx)'));

% Add a title...
title(titText);
    
box on;

set(gca,...
    'XTickLabel',[-fcThresh fcThresh],...
    'YTickLabel',{['q = ' num2str(10^-pqThresh)],...
        ['q = ' num2str(10^-floor(max(pq)))]},...
    'XTick',[-fcThresh fcThresh],...
    'YTick',[pqThresh floor(max(pq))],...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');
grid off;

% Draw lines along the origin
absFC = max(abs(fc)) * 1.1;
if absFC <= fcThresh
    xlim([-fcThresh-0.2 fcThresh+0.2]);
else
    xlim([-absFC absFC]);
end
    

%xlim('auto');
ylim([-0.1 max(pq) * 1.05]);

% Determine the current limits for the lines
xl = xlim(gca);
yl = ylim(gca);

% Horizontal line for PQ significance
line([xl(1) xl(2)],[pqThresh pqThresh],...
    'LineStyle','--',...
    'Color','k');

% Add in another for a higher value...
line([xl(1) xl(2)],[floor(max(pq)) floor(max(pq))],...
    'LineStyle',':',...
    'Color','k');

% Vertical lines for FC significance
line([-fcThresh -fcThresh],[yl(1) yl(2)],...
    'LineStyle','--',...
    'Color','k');
line([fcThresh fcThresh],[yl(1) yl(2)],...
    'LineStyle','--',...
    'Color','k');

% Display the significant ones
ff = array2table([mz(fx)' pq(fx) fc(fx)],...
    'VariableNames',{'m_z','neg_log10_q','log2_FC'});
disp(ff);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doStemPlot(parent,parent2,mz,fc,pq,titText,histID,sp,anno)
% SCatter plot the MVA results

% Significant thresholds
pqThresh = -log10(0.01);

% Log the pq values and use these throughout.  Need to remove the Inf
% values, or set them to be above the others by a slight way
lowPQ = min(pq(pq > 0));
idxLo = pq == 0;
pq(isnan(pq)) = 1;
pq = -log10(pq);
pq(idxLo) = -log10(lowPQ) * 1.2;

% Do below but for the plot above
f0 = get(parent2,'Children');
title(parent2,'');
delete(f0);

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);

% Hold the axes
axes(parent);
hold on;

% Highlight only significant ones
%fx = abs(fc) > fcThresh & pq > pqThresh;
fx = pq > pqThresh;

if sum(fx) == 0
    fx = false(size(pq));
end

% Scatter away - bad variables
scatter(mz(~fx),pq(~fx),60,[0.2627 0.5765 0.7647],'o','filled',...
    'MarkerEdgeColor',[0.7 0.7 0.7]);
    
% These are the interesting ones
scatter(mz(fx),pq(fx),120,[0.8392 0.3765 0.3020],'o','filled',...
    'MarkerEdgeColor',[0.7 0.7 0.7],...
    'HitTest','on',...
    'ButtonDownFcn',{@volcanoClick,...
    [parent parent2],mz(fx)',pq(fx),sp(:,fx),histID,mz(fx),anno(:,[3 5])});

% Add a title...
title(titText);
    
box on;

set(gca,...
    'YTickLabel',{['q = ' num2str(10^-pqThresh)],...
        ['q = ' num2str(10^-floor(max(pq)))]},...
    'YTick',[pqThresh max([floor(max(pq)) pqThresh+1])],...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal',...
    'XTickMode','auto',...
    'XTickLabelMode','auto');
grid off;

% Draw lines along the origin
xlim([min(mz) max(mz)]);
    
%xlim('auto');
ylim([-0.1 max(pq) * 1.05]);

% Determine the current limits for the lines
xl = xlim(gca);

% Horizontal line for PQ significance
line([xl(1) xl(2)],[pqThresh pqThresh],...
    'LineStyle','--',...
    'Color','k');

% Add in another for a higher value...
line([xl(1) xl(2)],[floor(max(pq)) floor(max(pq))],...
    'LineStyle',':',...
    'Color','k');

% Display the significant ones
try
    aa2 = [mz(fx)' pq(fx) fc(fx)];
catch
    aa2 = [mz(fx) pq(fx) fc(fx)];
end    
    
ff = array2table(aa2,...
    'VariableNames',{'m_z','neg_log10_q','log2_FC'});
disp(ff);
assignin('base','sigDiff',ff);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doLoadingsPlot(parent,mz,loadings,titText,histIDs)
% PLot the loadings

f0 = get(parent,'Children');
delete(f0);

axes(parent);
hold on;

% Insert zeros...
[x1,y1] = insertZeros(mz,loadings(:,1:2)',0.01);

% Apply an offset to make these visible
y1(2,:) = y1(2,:) - (0.5 * max(y1(2,:))); 

%plot(mz,loadings(:,1:2));
plot(x1,y1);

title(titText);

set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');

grid off;

xlim([min(x1) max(x1)]);
ylim([min(y1(:)) max(y1(:))]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doLoadingsPlotDual(parent,mz,loadings,fuseMethod,histIDs)
% PLot the loadings

f0 = get(parent,'Children');
delete(f0);

axes(parent);
hold on;

% Define the split between the values by looking for a drop in m/z values
fx = mz < [0 mz(1:end-1)];

% These are the two chunks of data
i1 = (1:numel(mz)) <  find(fx);
i2 = (1:numel(mz)) >= find(fx);

% Insert zeros... into the two separate parts...
[x1,y1] = insertZeros(mz(i1),loadings(i1,1:2)',0.01);
[x2,y2] = insertZeros(mz(i2),loadings(i2,1:2)',0.01);

% Apply an mz offset if showing high-level loadings, i.e. of scores
if strcmpi(fuseMethod(1:2),'hl')
    x1 = x1 - 0.05;
    x2 = x2 + 0.05;
end

% Apply an offset to make these visible
os = 0.5 * max([max(y1(:)) max(y2(:))]);
y1(2,:) = y1(2,:) - os; 
y2(2,:) = y2(2,:) - os;

% Plot the loadings based on those with more variables first
if numel(x1) > numel(x2)    
    plot(x1,y1,'Color',[0.5586 0.2188 0.7578]);
    plot(x2,y2,'Color',[0.4180 0.7578 0.2188]);    
else
    plot(x2,y2,'Color',[107 194 56]/256);
    plot(x1,y1,'Color',[143 56 194]/256);    
end

titText = ['\fontsize{11}Loadings of PC1 (top) and PC2 (bottom). '...
    '{\color[rgb]{0.5586 0.2188 0.7578}Positive}/'...
    '{\color[rgb]{0.4180 0.7578 0.2188}negative}.'];

title(titText);

set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');

grid off;

% Min/max y values?
minY = [min(y1(:)) min(y2(:))];
maxY = [max(y1(:)) max(y2(:))];

if strcmpi(fuseMethod(1:2),'hl')
    xlim([min(mz)-0.5 max(mz)+0.5]);
else
    xlim([min(x1) max(x1)]);
end

ylim([min(minY) max(maxY)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doSpectralPlot(parent,mz,sp,titText,histIDs,anno)
% Plot the spectra

f0 = get(parent,'Children');
delete(f0);

axes(parent);
hold on;

% Insert zeros...
[x1,y1] = insertZeros(mz,sp,0.01);
%y1 = sp;
%x1 = mz;

% Determine unique groupings...
[unq,~,ind] = unique(histIDs);
numG = numel(unq);

% Loop through
for n = 1:numG
    
    fx = ind == n;
        
    % Colour finding/matching
    cx = strcmp(anno(:,5),unq{n});
    cx = find(cx == 1,1,'first');
    col = anno{cx,3};

    if numG == 2 && n == 2
        mod = -1;
    else
        mod = 1;
    end
    
    plot(x1,y1(fx,:) * mod,'Color',col,...
        'LineWidth',2);
    
    
end

title(titText);

set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal');

grid off;

xlim([min(x1) max(x1)]);
if numG == 2
    nn = y1(fx,:);
    pp = y1(~fx,:);
    ylim([max(nn(:))*-1 max(pp(:))]);
else
    ylim([min(y1(:)) max(y1(:))]);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doHeatMap(parent,mz,sp,titText,histIDs,anno)
% Plot the spectra

% Sort the spectra according to histID
[~,srt] = sort(histIDs);
sp = sp(srt,:);
histIDs = histIDs(srt,:);

% Convert the spectra to be log2 fold changes
sp = bsxfun(@rdivide,sp,nanmedian(sp,1));
sp = log2(sp);
sp = real(sp);

sp(isnan(sp)) = 0;
sp(isinf(sp)) = 0;

f0 = get(parent,'Children');
delete(f0);

imagesc(sp,'Parent',parent);
hold on;

colormap(redbluecmap);

% What about the colour limits?
maxFC = max(abs(sp(:)));
caxis(parent,[-maxFC maxFC]);


% Draw some dividing lines...
[unq,vals,ind] = unique(histIDs);
numG = numel(unq);
for n = 1:numG
    
    % Colour finding/matching
    cx = strcmp(anno(:,5),unq{n});
    cx = find(cx == 1,1,'first');
    col = anno{cx,3};

    % Two places to draw the lines...
    st = vals(n);
    if n == numG
        fn = size(sp,1)+1;
    else
        fn = vals(n+1);
    end
    
    line([0 size(sp,2)],[st st]-0.5,...
        'LineWidth',3,...
        'Parent',parent,...
        'Color',col,...
        'LineStyle','--');
    
    line([0 size(sp,2)],[fn fn]-0.5,...
        'LineWidth',3,...
        'Parent',parent,...
        'Color',col,...
        'LineStyle','-');

    % Add a scatter point
    scatter(parent,1.5,mean([st fn]-0.5),200,col,'d','filled',...
        'MarkerEdgeColor','black');
    
    
end
    

% Axes properties
set(parent,'YTickLabel',[],...
    'YTick',[],...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal',...
    'XTick',[],...
    'Box','on',...
    'XLim',[0.5 size(sp,2)],...
    'YLim',[0.25 size(sp,1)+0.75]);
title(parent,'Heat map of fold changes relative to the overall median');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function volcanoClick(src,event,fig,lfc,qval,sp,histID,mz,cols)
% Get figure click... and then...

% Get coordinates
coor = get(fig(1),'CurrentPoint');
x = coor(1,1);
y = coor(1,2);

% Delete existing hitPoints...
f0 = findobj('Tag','hitPoint');
delete(f0);

%lfc(isinf(lfc)) = 0;

% Calculate the nearest to the lfc
dFC = (lfc - x) .^ 2; 
dFC = dFC / nanmax(dFC);
dQV = (qval - y) .^ 2; 
dQV = dQV / nanmax(dQV);
[~,dd] = min(dFC + dQV);

% Add a marker for the clicked point
hitPoint = scatter(lfc(dd),qval(dd),...
    600,[0.8392 0.3765 0.3020],'p','filled',...
    'MarkerEdgeColor',[0.8392 0.3765 0.3020]);
set(hitPoint,'Tag','hitPoint');

% Reset the secondary axes
axes(fig(2));
f0 = get(fig(2),'Children');
delete(f0);
%cla reset;

% Determine the unique values
[unq,~,ind] = unique(histID);

% Here we must determine the colours of the groups defined above
cidx = zeros(numel(unq),1);
for n = 1:numel(unq)    
    fx = strcmp(cols(:,2),unq{n});
    cidx(n,1) = find(fx,1,'first');    
end
allcols = vertcat(cols{:,1});
cols = allcols(cidx,:);
    

bpmethod = 'jsm';
switch bpmethod
    
    case 'jsm'

        % Use my default boxplot function
        jsmBoxPlot(sp(:,dd),ind,...
            'Orientation',fig(2),...
            'Colours',cols,...
            'Legend',false,...
            'Labels',unq,...
            'Order',1:numel(unq));
       
        % Change some of the axes properties
        set(fig(2),...
            'YTickLabel',[],...
            'FontSize',10,...
            'FontWeight','normal',...
            'FontName','Helvetica',...
            'YDir','normal',...
            'YTickMode','auto',...
            'YTickLabelMode','auto');
        
        ylabel(fig(2),'');
        %title(fig(2),'');

    case 'trad'

        % Now do the boxplot of this variable...
        boxplot(fig(2),...
            sp(:,dd),...
            ind,... % groups
            'Labels',unq);

end

txt = ['m/z = ' sprintf('%0.4f',mz(dd)) ...
    ', log2FC = ' sprintf('%0.2f',lfc(dd)) ...
    ', q = ' sprintf('%0.4e',10^-qval(dd))];

title(fig(2),txt,'FontSize',10);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
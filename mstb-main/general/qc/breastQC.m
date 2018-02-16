function [new] = breastQC(orig,fNames)
%breastQC - starting with the aligned & peak picked files, we need to get
%rid of files that look like outliers to leave us with good quality spectra

% Makes a lot of plots...
verbose = true;

% Extract the make names shorter
meta = orig.meta;
raw = orig.sp;

% Determine what class information we will use below to determine the
% groupings.
clInfo = classMany2One([meta.hist meta.cutcoag]);
%clInfo = meta.hist;

% First we will change the m/z range and re-normalise to account for the
% overbearing presence of the peak at 205.1
mask = orig.mz >= 0;
new.mz = orig.mz(mask);
new.sp = orig.sp(:,mask);

% Normalise
[new.sp] = jsmNormalise(new.sp,'pqn-median',0,0);

% Rescale to a similar level...
scFac = median(nansum(orig.sp(:,mask),2));
new.sp = new.sp * scFac;

% Log the data
logOS = nanmedian(new.sp(new.sp > 0));
unlog = new.sp;
new.sp = log(new.sp + logOS);

% Create a vector that says if we are keeping the files...
include = false(size(new.sp,1),1);

% Now we could proceed on a class by class basis to see what the
% differences are...
[unqCl,~,unqInd] = unique(clInfo);
numG = numel(unqCl);

if verbose
    figure;
    ax = scatter(1,1);
end

% Loop through each of the classes
for n = 1:numG

    % Observations in this group
    fx = unqInd == n;
    
    % Now many observations?
    fy = find(fx == 1);
    numO = numel(fy);
    
    % Can't do it for so few samples... so must include them
    if numO < 3
        include(fx,1) = true;
        continue;
    end        
    
    % Temporary data
    tmpD = new.sp(fy,:);
    
    % All mahal distances
    mhl = zeros(numO,1);
    
    pred = zeros(numO,size(tmpD,1)-1);
        
    % Loop though each observation
    for r = 1:numO
        
        hold off;
        
        % Determine which to be included
        inc = fy ~= fy(r);
        
        % Run PCA
        [ll,ss,~,~] = princomp(tmpD(inc,:),'econ');
        
        % Which components to use?
        comps = [1 2]; %[size(ss,2)-1 size(ss,2)];
                
        % Determine 95% confidence ellipse...
        [ellPl,centr,ellAx] = error_ellipse(ss(:,comps),95);
                
        % Predict new sample
        pred(r,1:size(ll,2)) = bsxfun(@minus,tmpD(~inc,:),nanmean(tmpD(inc,:),1)) * ll;
        
        % Calculate the Mahalanobis distance
        mhl(r,1) = mahal(pred(r,1:size(ll,2)),ss);
        
        % Check point
        c1 = (pred(r,comps(1))-centr(1)) ^ 2 / ellAx(1)^2;
        c2 = (pred(r,comps(2))-centr(2)) ^ 2 / ellAx(2)^2;
        
        % This is the inclusion / exclusion criteria.
        if c1 + c2 <= 1
            flag = true;
            include(fy(r),1) = true;
        else
            flag = false;
            include(fy(r),1) = false;
        end

        if verbose && ~flag
            cla reset
            hold on;
            scatter(ss(:,comps(1)),ss(:,comps(2)),80,'k','filled');
            scatter(pred(r,comps(1)),pred(r,comps(2)),200,'r','filled');
            plot(ellPl(:,1),ellPl(:,2));
            drawnow;
            pause(0.01);
        end
                
    end
    
    % Now with the list of excluded samples from this group, we need to
    % show how different they are from the rest of the samples
    keep = include(fy);
        
    % So now scatter plot the predicted PCA scores from the leave one out
    % approach, and colour the excluded samples accordingly.  Put an
    % ellipse round in order to show the outliers.
    figure; hold on;
    scatter(pred(keep,comps(1)),pred(keep,comps(2)),80,'b','o','filled');
    scatter(pred(~keep,comps(1)),pred(~keep,comps(2)),80,'r','o','filled');
    [ellPl,~,~] = error_ellipse(pred(keep,comps),95);
    plot(ellPl(:,1),ellPl(:,2),'LineWidth',2);
    box on;    
    
    set(gca,'FontSize',14);
    xlabel('PC1 - Projection','FontSize',16,'FontWeight','bold');
    ylabel('PC2 - Projection','FontSize',16,'FontWeight','bold');
    title(unqCl{n},'FontSize',16,'FontWeight','bold')
    
    graphFormat([pwd filesep unqCl{n} '-PCA'],'png');
    close
    
    % What about spectral profiles...
    if sum(keep) == numel(keep)
        continue;
    end
    
    [t1,t2] = insertZeros(new.mz,nanmean(unlog(fx & include,:),1));
    [u1,u2] = insertZeros(new.mz,unlog(fx & ~include,:));
    figure; hold on;
    plot(u1,u2,'r','LineWidth',2);
    plot(t1,t2,'b','LineWidth',1);
    graphFormat([pwd filesep unqCl{n} '-Plot'],'png');
    close
    
end

% Unlog the data so that it is back to normal...
new.sp = unlog;

% Keep the 'old' meta data
new.oldMeta = meta;

% Create a new metadata variable saying which were included/excluded...
meta.outlier = cell(size(include));
meta.outlier(include,1) = {'Normal'};
meta.outlier(~include,1) = {'Outlier'};

% Trim the dataset and metadata to exclude the desired observations, but do
% it from the raw data so that normalisation can be performed again
% differently later on as desired
[new.sp,new.meta] = removeObservations(raw,meta,~include);
new.mz = orig.mz;

% List the vector of files that were included as it is useful to know
new.include = include;
new.files = fNames(include);

% THis is how the classes were divided
new.classes = clInfo;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quickPCA(x,y)
% Do PCA and plot

[~,s,ee] = princomp(x,'econ');

[ ax,centr ] = scatEll([],s(:,1:2),y,[1 2])
[ ax,centr ] = scatEll([],s(:,end-1:end),y,[size(s,2)-1 size(s,2)])

%[ hFig,hAx ] = pca_ell(s,ones(size(y)),y);

return

scatterPlot(s(:,1:2),y,[],'PC1','PC2');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spectraPlot(mz,x,y)

[unq,~,ind] = unique(y);
cols = jet(numel(unq));

figure; hold on;
for n = 1:numel(unq)    
    fx = ind == n;    
    plot(mz,x(fx,:),'Color',cols(n,:));    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
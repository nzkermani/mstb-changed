function statsTSNEClusters(src,event,fig,man)
%statsTSNEClusters - calculate clusters based on the tsne data and the
%options provided...

% Guidata
sts = guidata(fig.fig);
if ~isfield(sts,'res')
    return
elseif ~isfield(sts.res,'tsne')
    return
end

% Get the clustering options
valDist = man.clustDist.Value;
if valDist == 1
    valDist = 2;
    man.clustDist.Value = 2;
end
opts.clustDist = man.clustDist.String{valDist};

valLink = man.clustLink.Value;
if valLink == 1
    valLink = 2;
    man.clustLink.Value = 2;
end
opts.clustLink = man.clustLink.String{valLink};

opts.clustMax = str2double(man.clustMax.String);

% Generate distance matrix
dst = pdist(sts.res.tsne.ss,opts.clustDist);

% Link the distances
lnk = linkage(dst,opts.clustLink);

% Cophenetic
cc = cophenet(lnk,dst);
disp(['Cophenetic correlation = ' sprintf('%0.2f',cc)]);

% Finalise clustering
t = cluster(lnk,'MaxClust',opts.clustMax);

% What is the histological group that we are interested in?
val = man.groups.Value;
str = man.groups.String;
grp = statsObservationLabels(sts.proc.meta,str,val);

% Now plot the results in the top right axes
[patchHand] = statsScatterPlotWithClusters(fig.ax.conf(1),sts.res.tsne.ss(:,1:2),grp,...
    false,false,sts.proc.meta,t);

% Now we need to save the clusters and the handldes to the patches in order
% that we can do something fun with them...
sts.res.tsne.clusts = t;
sts.res.tsne.patch = patchHand;

% Add callback functions to the patches
set(sts.res.tsne.patch,'ButtonDownFcn',{@statsTSNEpatchCallback,fig,man});

% I also want to write the cluster number into the 'raw' part of meta data
% so that we can go and remove those observations at a later time, or just
% focus on a specific cluster or two
sts.proc.meta.clustID = t;

% Match cluster assignments
[sts.raw] = matchObs2Clust(sts.raw,sts.proc);

% Then we need to update the table to include the clusters
statsTablePopulate([],[],fig);

% Save guidata
guidata(fig.fig,sts);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raw] = matchObs2Clust(raw,proc)
% Match any cluster labels obtained from this 'proc' version of the data
% and ensure that they are matched to the appropriate parts in the 'raw'
% meta data.  In instances where all data was used

% In instances where all observations survived from raw->proc, then no
% matching needs to be performed, in theory, as there is a direct
% comparison exists.  When obs(proc) < obs(raw), then we need to match

% Easiest to do in a loop?
rawClust = zeros(size(raw.obsID));

% What are the IDs of the observations in proc?
procID = proc.obsID;

% Add in...
rawClust(procID) = proc.meta.clustID;

% Save in raw
raw.meta.clustID = rawClust;

% We also need to update the table to show the clustIDs

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function statsTSNEpatchCallback(src,~,fig,man)
%statsTSNEpatchCallback - do stuff when people click a patch

% Get the guidata
sts = guidata(fig.fig);

% Which cluster did we click
thisC = str2double(src.Tag(6:end));
disp(['Clicked Cluster = ' int2str(thisC)]);

% Cluster labels from tsne
t = sts.res.tsne.clusts;

% Observations in these clusters
obs = t == thisC;

% Which group are we drawing as coloured in bits on the plots?
grp = man.groups.String{man.groups.Value};

% Simple analysis of observations in this cluster...
thisFold = sts.proc.meta.foldID(obs);

dispGroup = sts.proc.meta.(grp);
thisGroup = dispGroup(obs);

freqHist = cellFreq(thisGroup);
freqFold = cellFreq(thisFold);

end


function volcanoPlot2(mz,sp,grp,qVal,fldChange)
%volcanoPlot2 -  a standalone version of the original volcano plot

qVal = -log10(qVal);

fldChange(isinf(fldChange)) = 0;
fldChange(isnan(fldChange)) = 0;


% Do some basic housekeeping:
% - set all 0-fold change variables to qVal = NaN;
% - set all 0-qVal variables to fold change = NaN;
noFC = fldChange == 0;
noQV = qVal == 0;
fldChange(noQV) = 0;
qVal(noFC) = 0;

% Create a guidata structure to make the figure easier to use...
g.mz = mz;
g.sp = sp;
g.grp = grp;
g.groupIds = unique(g.grp);
g.qVal = qVal;
g.fldChange = fldChange;

% Draw a figure...
fig.fig = figure('Name','Volcano!'); hold on;
guidata(fig.fig,g);

% Draw the subaxes...
fig.ax(2) = subplot(1,3,2);
fig.ax(3) = subplot(1,3,3);

[mz2,sp2] = insertZeros(mz,sp,0.01);
plot(mz2,sp2);

fig.ax(1) = subplot(1,3,1);

% Draw the main volcano plot
scatter(fldChange,qVal,40,mz,'o','filled',...
    'MarkerEdgeColor','k',...
    'HitTest','on',...
    'ButtonDownFcn',{@figClick,fig,fldChange,qVal}); 

%cmap = flipud(gray(100));
%colormap(cmap);

title(['Fold changes'],...
    'FontSize',16);
xlabel('log_2 fold change','FontSize',16);
ylabel('-log_1_0 q-value','FontSize',16);
cbar = colorbar;
ylabel(cbar,'m/z','FontSize',16);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figClick(src,event,fig,lfc,qval)
% Get figure click... and then...

% Get coordinates
coor = get(fig.ax(1),'CurrentPoint');
x = coor(1,1);
y = coor(1,2);


% Get the actual data...
g = guidata(fig.fig);

% Calculate the nearest to the lfc
dFC = (lfc - x) .^ 2; 
dFC = dFC / max(dFC);
dQV = (qval - y) .^ 2; 
dQV = dQV / max(dQV);
[~,dd] = min(dFC + dQV);

% Reset the secondary axes
axes(fig.ax(2));
%cla reset;

% Now do the boxplot of this variable...
boxplot(fig.ax(2),g.sp(:,dd),g.grp);%,'Labels',g.groupIds);
ylabel(fig.ax(2),'Intensity','FontSize',16);

txt = ['m/z = ' sprintf('%0.4f',g.mz(dd)) ...
    ', log_2FC = ' sprintf('%0.2f',lfc(dd)) ...
    ', q = ' sprintf('%0.4e',10^-qval(dd))];
title(fig.ax(2),txt,'FontSize',16);

set(fig.ax(1),'ButtonDownFcn',{@figClick,fig,lfc,qval});%,'HitTest','off');


% Display statistics for each of the groups...
disp(sprintf('%0.4f',g.mz(dd)));
[unqG,~,indG] = unique(g.grp);
for n = 1:numel(unqG)
    
    disp(unqG{n});
    intVals = prctile(g.sp(indG == n,dd),[25 50 75]);
    disp(intVals);
end
disp('-----------------------------');


% Now draw the peaks in the third axes to give an idea if it is actually a
% peak or noise or something else
axes(fig.ax(3));
xlim([g.mz(dd)-1.1 g.mz(dd)+1.1]);
ylim([-1 max(g.sp(:,dd))*1.25]);

drawnow;




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
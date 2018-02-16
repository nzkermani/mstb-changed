function mtspDiagnostics(info)
%mtspDiagnostics - make some fancy plots which show how the annotation,
%recalibration and reannotation has performed.  Specifically about numbers
%of annotations...

% Which files were completed successfully?
st = {info.status}';
fx = strcmp(st,'Completed');

% List the files with no toolbox file
none = info(~fx);
disp('UNFINISHED FILES');
disp({none.file}');

% Now we can focus on the ones for which we have information
info = info(fx);

% Combine numbers
annos = [vertcat(info.qtyOld) vertcat(info.qtyNew) vertcat(info.msiOld) vertcat(info.msiNew)];

% File names
fnames = {info.file}';

% Deviations
devs = vertcat(info.precal);

% Draw a figure
annoPlotMTSP(annos(:,1),annos(:,2),median(devs,2),...
    'Annotations (original)',...
    'Annotations (recalibrated)',...
    'Median ppm deviation',...
    fnames,true);

annoPlotMTSP(annos(:,3),annos(:,4),median(devs,2),...
    'MSI features (original)',...
    'MSI features (recalibrated)',...
    'Median ppm deviation',...
    [],false);

annoPlotMTSP(annos(:,3),annos(:,1),median(devs,2),...
    'MSI features (original)',...
    'Annotations (original)',...
    'Median ppm deviation',...
    [],false);

annoPlotMTSP(annos(:,4),annos(:,2),'k',...
    'MSI features (recalibrated)',...
    'Annotations (recalibrated)',...
    'Median ppm deviation',...
    [],false);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function annoPlotMTSP(x,y,c,xlab,ylab,clab,labels,equal)

fx = x > y;

figure('Units','pixels','Position',[642 594 918 744]); hold on;


scatter(x,y,80,c,'o','filled',...
    'MarkerEdgeColor',[0.5 0.5 0.5]);

if ~isempty(labels)

    scatter(x(fx),y(fx),200,c(fx),'s','filled',...
        'MarkerEdgeColor',[0 0 0]);

    text(x(fx)+10,y(fx),labels(fx),'FontSize',14);
end

if equal
    axis equal
    xl = xlim;
    yl = ylim;
    mv = max([xl yl]);
    line([0 mv],[0 mv]);
    xlim([0 mv]);
    ylim([0 mv]);
end

box on;
set(gca,'FontSize',14);

xlabel(xlab,'FontWeight','bold','FontSize',16);
ylabel(ylab,'FontWeight','bold','FontSize',16);

if isnumeric(c)
    cmap = flipud(redgreencmap(101));

    cb = colorbar;
    colormap(cmap);
    ylabel(cb,clab,'FontSize',16,'FontWeight','bold');

    cl = max(abs(c));
    caxis([-cl cl]);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
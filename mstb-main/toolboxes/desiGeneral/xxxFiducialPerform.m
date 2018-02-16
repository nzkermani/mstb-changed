function xxxFiducialPerform(src,event,fig,man)
%xxxFiducialPerform - carry out the actual coregistration using the
%fiducial markers

% Get the guidata
dpn = guidata(fig.fig);
if ~isfield(dpn.opt,'fiducial')
    return
elseif isempty(dpn.opt.fiducial)
    return
end

% How many points are there?
sz = size(dpn.opt.fiducial,1);
if sz < 3
    disp('Must have at least 3 points for transformation');    
    return
end

% Which method are we using?
method = get(man.method,'String');
vals = get(man.method,'Value');
method = method{vals};

% Get the reference image
refImage = get(fig.ax.ms1(2),'CData');

% Need a reference size based on the master image
Rfixed = imref2d(size(refImage));
%RfixedL = imref2d(size(dpn.d1.sp) * 10)

% Now we need to determine which points are master (MS) and which points
% are to be moved (optical)
ptOP = cell2mat(dpn.opt.fiducial(:,7:8));
ptMS = cell2mat(dpn.opt.fiducial(:,9:10));

% Geometric transformation
transform = fitgeotrans(ptOP,ptMS,method);

% Now warp the optical image
visOptImg = get(fig.ax.opt(2),'CData');
[dpn.opt.coreg,~] = imwarp(visOptImg,transform,'OutputView',Rfixed);

% And transform the points as well
[nx,ny] = transformPointsForward(transform,ptOP(:,1),ptOP(:,2));


% Create the figure to show success / failure of the coregistration
makeFigure(dpn.opt.fiducial,ptMS,nx,ny,dpn.opt.coreg,refImage);

% Only continue if we want to finalise the coregistration - may just want
% to try out the various methods...
fin = get(man.finalise,'Value');
if ~fin
    disp('Not finalised');
    return
end


% Now need to finish off the various steps, i.e. update the optical image,
% delete all of the fiducial markers and allow things like annotation to
% proceed...

% Delete the points from the axes
delete([dpn.opt.fiducial{:,5}]);
delete([dpn.opt.fiducial{:,6}]);

% And actually delete the fiducial markers because we are done with them
dpn.opt.fiducial = [];

% Empty the table of points
set(man.tab,'Data',[]);

% Save this information to the structure...
guidata(fig.fig,dpn);

% Update the plots to show all the right stuff again...
set(fig.ax.opt(2),'CData',dpn.opt.coreg,...
    'XData',[1 size(dpn.opt.coreg,2)],...
    'YData',[1 size(dpn.opt.coreg,1)]);
set(fig.ax.opt(1),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.coreg,2)],...
    'YLim',[1 size(dpn.opt.coreg,1)]);

% Calculate the grid (again) to adapt to the (new) co-registration
desiDetermineGrid([],[],fig);

% Update both images...
dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);

% Then turn off the fiducial marker thing
set(fig.tb.fiducial,'State','off');
f0 = get(fig.pan2,'Children');
delete(f0);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeFigure(tab,ptMS,nx,ny,registered,refImage)
% Make a new figure that appears, which shows the success of the
% coregistration process...

% Let's determine the outline of the optical image...
ol = imageOutline(registered,false);

cols = vertcat(tab{:,4});

figure('Units','normalized','Position',[0.25 0.25 0.5 0.5]);

ax1 = subplot(1,2,1); hold on;
imagesc(registered);
scatter(nx,ny,200,cols,'o','filled',...
    'MarkerEdgeColor','k');
scatter(ptMS(:,1),ptMS(:,2),80,cols,'d','filled',...
    'MarkerEdgeColor','k');

plot(ol(:,2),ol(:,1),'LineWidth',2,'Color','k');

set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',[]);
xlim([0.5 size(registered,2)+0.5]);
ylim([0.5 size(registered,1)+0.5]);

ax2 = subplot(1,2,2); hold on;
imagesc(refImage);
colormap redbluecmap
scatter(ptMS(:,1),ptMS(:,2),200,cols,'d','filled',...
    'MarkerEdgeColor','k');
scatter(nx,ny,80,cols,'o','filled',...
    'MarkerEdgeColor','k');

plot(ol(:,2),ol(:,1),'LineWidth',2,'Color','k');

set(gca,'YDir','reverse','XTickLabel',[],'YTickLabel',[]);
xlim([0.5 size(registered,2)+0.5]);
ylim([0.5 size(registered,1)+0.5]);

% Link the axes
linkaxes([ax1 ax2],'xy');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xxxFiducialReset(~,~,fig,man)
%xxxFiducialReset - restore the original optical image, remove all fiducial
%markers, and prepare to start all over again...

% Get guidata
dpn = guidata(fig.fig);

% Get markers
if isfield(dpn.opt,'fiducial')
    fid = dpn.opt.fiducial;
end
% if isempty(fid)
%     return
% end

% Delete all scatter points
if ~isempty(fid)
    delete([fid{:,5}]);
    delete([fid{:,6}]);
end

% Reload optical image
set(fig.ax.opt(2),'CData',dpn.opt.orig);
set(fig.ax.opt(1),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dpn.opt.orig,2)],...
    'YLim',[1 size(dpn.opt.orig,1)]);

% Delete markers from dpn
dpn.opt.fiducial = [];

% Remove all entries from the table
set(man.tab,'Data',[]);

% Update guidata
guidata(fig.fig,dpn);

% Return the upsized MS image
tmpImg = imresize(dpn.d1.img,10);
dpnIonImage([],[],fig.ax.ms1,tmpImg);


end


function dpnCallbacks(~,~,fig,defP)
%toolbarCallbacks - update the toolbar's callback functions in here

% This is for reloading an already analysed dataset
set(fig.tb.new,'ClickedCallback',{@dpnIMZML,fig,defP});

% Set the callback for the change path button
set(fig.tb.save,'ClickedCallback',{@dpnSave,fig});

% Set the callback for the go button
set(fig.tb.refresh,'ClickedCallback',{@dpnUpdateMS,fig});

% Set the callback to switch the images
set(fig.tb.switch,'ClickedCallback',{@dpnSwitch,fig});

% Set the callback for the go button
set(fig.tb.optimg,'ClickedCallback',{@dpnOptImg,fig});

% Set the pos/neg manipulations
set(fig.tb.editPos,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@dpnManipulateMS,fig,'pos'}});
set(fig.tb.editNeg,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@dpnManipulateMS,fig,'neg'}});

% Image interpolation
% set(fig.tb.interp,'ClickedCallback',...
%     {@desiDrawSideMenu,fig,{@dpnInterpolateCallback,fig}});
set(fig.tb.interp,'ClickedCallback',{@dpnInterpolateCallback,fig});

% Image coregistration
set(fig.tb.coreg,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@dpnCoreg,fig}});

% Zoom
set(fig.tb.zoom(1),'ClickedCallback',{@desiZoomCallback,fig,+1});
set(fig.tb.zoom(2),'ClickedCallback',{@desiZoomCallback,fig,-1});
set(fig.tb.zoomreset,'ClickedCallback',{@desiZoomReset,fig});

% Grid
set(fig.tb.grid,'ClickedCallback',{@desiGrid,fig});

% Annotation
set(fig.tb.annotate,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@dpnAnnotate,fig}});

% Internal statistics
set(fig.tb.intStats,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@desiStatistics,fig}});


% Segmentation
set(fig.tb.segment,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@dpnSegmentationWindow,fig}});

% Segmentation
set(fig.tb.family,'ClickedCallback',...
    {@desiDrawSideMenu,fig,{@xxxLipidFamilyWindow,fig}});

% Fusion
%set(fig.tb.fuse,'ClickedCallback',{@dpnFusionCallback,fig});

% Layout
set(fig.tb.layout,'ClickedCallback',{@dpnChangeLayout,fig});

% Correlation comparison
set(fig.tb.doComparison,'ClickedCallback',{@dpnComparisonCallback,fig});

% Export figures
set(fig.tb.expfig,'ClickedCallback',{@xxxExportFigures,fig});

% Close all other windows
set(fig.tb.duck,'ClickedCallback',@closeOtherWindows);


end


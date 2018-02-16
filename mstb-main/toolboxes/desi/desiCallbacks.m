function desiCallbacks(~,~,fig,defP)
%desiCallbacks - update the toolbar's callback functions in here

% This is for reloading an already analysed dataset
%set(fig.tb.new,'ClickedCallback',{@desiFileUpload,fig,defP});
set(fig.tb.new,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@desiFileUpload,fig,defP}});

% Triple Q import
set(fig.tb.qqq,'ClickedCallback',{@desiQQQUpload,fig,defP});

% Set the callback for the change path button
set(fig.tb.save,'ClickedCallback',{@dpnSave,fig});

% Set the callback for the go button
set(fig.tb.refresh,'ClickedCallback',{@dpnUpdateMS,fig});

% Set the callback for the crop button
set(fig.tb.crop,'ClickedCallback',{@desiCrop,fig});

% Set the callback for the crop button
set(fig.tb.reproc,'ClickedCallback',{@desiReproc,fig});

% Show the QC images and decide if it is pass or fail
set(fig.tb.qcImage,'ClickedCallback',{@desiShowQC,fig});

% Set the callback for the crop button
%set(fig.tb.raman,'ClickedCallback',{@desiRaman,fig});

% Set the callback for polymer removal
%set(fig.tb.removePoly,'ClickedCallback',{@desiRemovePolymer,fig});

% Set the callback for the flipping
%set(fig.tb.flipUD,'ClickedCallback',{@xxxFlipUpDown,fig});

% Set the callback for the go button
set(fig.tb.optimg,'ClickedCallback',{@dpnOptImg,fig});

% Set the pos/neg manipulations
set(fig.tb.editPos,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@dpnManipulateMS,fig,'pos'}});

% Image coregistration
set(fig.tb.coreg,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@dpnCoreg,fig}});

% Fiducial coregistration
set(fig.tb.fiducial,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@xxxFiducial,fig}});

% Erasing lines in an image
set(fig.tb.erase,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@xxxErase,fig}});

% Enable histo image as an overlay
set(fig.tb.hevis,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@xxxHistoVis,fig}});

% Set the callback for flipping
% set(fig.tb.flipud,'ClickedCallback',{@xxxFlip,fig,'ud'});
% set(fig.tb.fliplr,'ClickedCallback',{@xxxFlip,fig,'lr'});


% Zoom in and out
set(fig.tb.zoom(1),'ClickedCallback',{@desiZoomCallback,fig,+1});
set(fig.tb.zoom(2),'ClickedCallback',{@desiZoomCallback,fig,-1});
set(fig.tb.zoomreset,'ClickedCallback',{@desiZoomReset,fig});

% Grid
set(fig.tb.grid,'ClickedCallback',{@desiGrid,fig});

% Annotation
set(fig.tb.annotate,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@dpnAnnotate,fig}});

% Internal statistics
set(fig.tb.intStats,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@desiStatistics,fig}});

% Stats or similar
%set(fig.tb.extract,'ClickedCallback',{@dpnDataExtract,fig});

% Segmentation
set(fig.tb.segment,...
    'ClickedCallback',{@desiDrawSideMenu,fig,{@dpnSegmentationWindow,fig}});

% Outlines
% set(fig.tb.outline,...
%     'ClickedCallback',{@desiDrawSideMenu,fig,{@desiOutlineWindow,fig}});

% Help
set(fig.tb.help,'ClickedCallback',{@desiHelp});


% Layout
set(fig.tb.layout,'ClickedCallback',{@desiChangeLayout,fig});

% Export figures
set(fig.tb.expfig,'ClickedCallback',{@xxxExportFigures,fig});

% Close all other windows
set(fig.tb.duck,'ClickedCallback',@closeOtherWindows);

end


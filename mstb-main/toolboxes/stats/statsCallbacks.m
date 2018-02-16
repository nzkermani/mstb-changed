function statsCallbacks(~,~,fig,defP)
%statsCallbacks - update the toolbar's callback functions in here

% Import from saved file
set(fig.tb.open,'ClickedCallback',{@statsOpen,fig,defP});

% Import instantly from workspace
set(fig.tb.import,'ClickedCallback',{@statsWorkspace,fig});

% Save the variables
set(fig.tb.save,'ClickedCallback',{@statsSave,fig,defP});

% Toggle the table (or not)
set(fig.tb.table,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsLayoutChange,fig}});

% Norm/tran/scal
set(fig.tb.norm,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsNormMenu,fig}});

% Multivariate
% set(fig.tb.statsMV,...
%     'ClickedCallback',{@statsSideMenuDraw,fig,{@statsMultivariateMenu,fig}});

% PCA
set(fig.tb.doPCA,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsPCAMenu,fig}});

% Cluster analyses
set(fig.tb.doClust,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsClusterMenu,fig}});

% tSNE
set(fig.tb.doTSNE,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsTSNEMenu,fig}});

% MMC/LDA/etc
set(fig.tb.doLDA,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsLDAMenu,fig}});

% External samples
set(fig.tb.doExternal,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsExternalMenu,fig}});

% Univariate
set(fig.tb.doUnivariate,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsUnivariateMenu,fig}});

% Univariate
set(fig.tb.doBoxplots,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsBoxplotsMenu,fig}});

% Ratios
set(fig.tb.doRatio,...
    'ClickedCallback',{@statsSideMenuDraw,fig,{@statsRatioMenu,fig}});

return

% Export figures
set(fig.tb.expfig,'ClickedCallback',{@xxxExportFigures,fig});

% Close all other windows
set(fig.tb.duck,'ClickedCallback',@closeOtherWindows);

end


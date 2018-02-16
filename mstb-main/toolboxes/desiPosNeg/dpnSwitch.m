function dpnSwitch(~,~,fig)
%dpnSwitch - change the pos/neg ion images around...

dpn = guidata(fig.fig);

% Determine the old ones
old1 = dpn.d1;
old2 = dpn.d2;

% Switch the things round
dpn.d1 = old2;
dpn.d2 = old1;

% Update the guidata
guidata(fig.fig,dpn);

% Update the images...
dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);
dpnIonImage([],[],fig.ax.ms2,dpn.d2.img);

end


function xxxFlipUpDown(src,event,fig)
%xxxFlipUpDown - mirror the DESI data

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return;
end

% Flip it
dpn.d1.sp = flipud(dpn.d1.sp);

% Save it
guidata(fig.fig,dpn);

% Update the ion image
dpnIonImage([],[],fig.ax.ms1,dpn.d1.img);

end


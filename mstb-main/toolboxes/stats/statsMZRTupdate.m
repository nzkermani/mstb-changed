function statsMZRTupdate(~,~,fig,man)
%statsMZRTupdate - change the entries in the window

sts = guidata(fig.fig);

val = man.mzrtchoice.Value;

op = statsMZRTRTMZ(sts.raw.var.mz,sts.raw.var.rt,val);

set(man.list,'String',op,'Value',1);

end


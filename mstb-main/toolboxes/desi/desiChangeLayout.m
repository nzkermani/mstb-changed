function desiChangeLayout(src,~,fig)
%desiChangeLayout - switch between the two layout modes

% Get the state of the button
value = get(src,'State');

switch value
    
    case 'on' % 'detail'
        % ...if on the change to the max detail view
        [locn] = desiLocations('detail');       
        
    case 'off' % 'simple'
        % ...if off change to the two-grid view
        [locn] = desiLocations('simple');
        
end

% Now run through the relevant axes and change the positions...
set(fig.ax.opt(1),'Position',locn.opt.ax);
set(fig.ax.ms1(1),'Position',locn.ms1.ax);

% Top right plot
set(fig.ax.sp(1), 'Position',locn.sp.ax);
set(fig.ax.sp,'Visible',locn.sp.vis,'XTick',[],'YTick',[]);
set(get(fig.ax.sp,'Children'),'Visible','off');
title(fig.ax.sp,'');


% Bottom right
set(fig.ax.mv(1), 'Position',locn.mv.ax);
set(fig.ax.mv,'Visible',locn.mv.vis,'XTick',[],'YTick',[]);
set(get(fig.ax.mv,'Children'),'Visible','off');
title(fig.ax.mv,'');

set(fig.titOpt,'Position',locn.opt.lab1,'Visible','off');
set(fig.titPos,'Position',locn.ms1.lab1,'Visible','off');
set(fig.txtOpt,'Position',locn.opt.lab2,'Visible','off');
set(fig.txtPos,'Position',locn.ms1.lab2,'Visible','off');



end


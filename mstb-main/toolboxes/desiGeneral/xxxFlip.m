function xxxFlip(src,event,fig,direction)
%xxxFlip - flip an image either LR or UD

% Gui and quit
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Simple definition of which flip to do
switch direction
    case 'ud'
        funcHand = @flipud;
        
    case 'lr'
        funcHand = @fliplr;        
end

dpn.d1.sp = funcHand(dpn.d1.sp);
dpn.d1.sp = funcHand(dpn.d1.sp);
dpn.d1.sp = funcHand(dpn.d1.sp);




end


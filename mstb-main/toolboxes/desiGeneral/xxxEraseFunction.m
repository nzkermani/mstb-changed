function xxxEraseFunction(src,event,fig,erase)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Get the number of lines and then the replacement value
lines = str2double(get(erase.lines,'String'));
value = str2double(get(erase.value,'String'));

% Set values
dpn.d1.sp(:,1:lines,:) = value;

% Update the guidata
guidata(fig.fig,dpn);

end


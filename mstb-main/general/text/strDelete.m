function [ new ] = strDelete( txt,ch )
%strDelete - delete all spaces from a string

%fx = isstrprop(txt,'wspace');
%fx = strfind(txt,ch);
fx = txt == ch;

new = txt(~fx);

end


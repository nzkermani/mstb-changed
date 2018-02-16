function [ new ] = strBlank( txt )
%strBlank - remove all spaces from a string

fx = isstrprop(txt,'wspace');

new = txt(~fx);

end


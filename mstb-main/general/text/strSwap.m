function [ txt ] = strSwap(txt,a,b)
%strSwap - replace 'a' in txt with 'b'

chk = strfind(txt,a);

txt(chk) = b;



end


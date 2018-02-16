function [elem] = parseFormula(form)
% Convert a string of chemical elements into a single homogeneous string
% that you can use in the function. Will require a defined list of
% elements, in order...
%
% James McKenzie, 2014.

% Define the elements
[elem] = defineElements;

% Shove a 1 on the end if absent...
if isstrprop(form(end),'lower') || isstrprop(form(end),'upper')
    form(end+1) = '1';
end

% Calculate the population: text / numeric
txt = isstrprop(form,'upper');
low = isstrprop(form,'lower');
num = isstrprop(form,'digit');


% How many chars in the formula?
lng = length(form);

% How many element entries are there?
numE = sum(txt);
ex = find(txt == 1);

for n = 1:numE
    
    % This is the element
    st = ex(n);    
    if low(ex(n)+1) == 1
        fn = ex(n)+1;
    else
        fn = st;
    end
    ee = form(st:fn);
    
    % Now the numbers... What is beyond fn?
    if num(fn+1) == 0
        qq = 1;
    else
        jj = fn+1;
        while num(jj) == 1 && jj < lng %|| jj <= lng-1
            jj = jj + 1;
        end
        
        if jj ~= lng
            jj = jj - 1;
        end
        qq = str2double(form(fn+1:jj));
    end
    
    % Now need to put into the qty vector in elem.qty
    ei = strcmp(elem.symb,ee);
    if sum(ei) == 0
        disp([ee ' is not an element that I recognise']);
    end
    elem.qty(ei) = elem.qty(ei) + qq;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

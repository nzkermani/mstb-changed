function [ mz ] = mtspForm2MZ(annos)
%mtspForm2MZ - convert annotations with adducts into m/z values

% Can provide the adduct as either part of the cell or separately...


numA = size(annos,1);

mz = zeros(numA,1);

for n = 1:numA
    
    % Parse formula
    elem = parseFormula(annos{n,1});
    
    % What is the adduct?
    switch annos{n,2}
        case '-H'
            elem.qty(2) = elem.qty(2) - 1;
        case '+Cl'
            elem.qty(12) = elem.qty(12) + 1;
    end
    
    % Now determine m/z value
    mz(n,1) = sum(elem.qty .* elem.mono);
    
end


end


function [ nom,kmd ] = kendrick(mz,sp,unit)
%kendrick - convert a series of m/z values into Kendrick values. Note that
%we will probably want to convert the unit from CH2 to something else

% Which repeat unit do we want to use?
switch lower(unit)
    case 'ch2'
        rto = 14 / 14.01565;
    case 'ch'
        rto = 13 / 13.007825;
    case {'c2h4o','peg'}
        rto = 44 / 44.026215;
    otherwise
        error('No other repeat units');
end

% Determine nominal unit integer mass (floor)
nom = floor(mz);

% Defect
kmd = nom - (mz * rto);

return

% Now plot it...
figure;
scatter(nom,kmd,80,sp,'o','filled');

figure;
stem3(nom,kmd,sp);



end


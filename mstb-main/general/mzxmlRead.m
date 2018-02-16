function mzxmlRead(fn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(fn)
    fn = '/Volumes/JSM/DB/Dried Blood Spots/MzZML Files by Patient/Cancer/Pt 39/TW03_H03.mzXML';
end

spec = mzxmlread(fn,'verbose','F');
numS = size(spec.scan,1);
tmp = cell(numS,3);

figure; hold on;



% Convert into something more resembling a series of mass spectra
for r = 1:numS
    
    % Get the retention time
    rt = spec.scan(r).retentionTime;
    fx = isstrprop(rt,'digit');
    fx = find(fx == 1);
    tmp{r,1} = str2double(rt(fx(1):fx(end)));
    
    % m/z and intensity vectors
    mz = spec.scan(r).peaks.mz(1:2:end);
    sp = spec.scan(r).peaks.mz(2:2:end);
    
    stem(mz,sp,'MarkerSize',0.1);
    
end



end


function [ output_args ] = pdmLocate( input_args )
%pdmLocate - find files in one folder and move to another

raw = '/Users/jmckenzi/Documents/Box Sync/CRUK GC/Colorectal Xevo/Negative/';
pdm = '/Users/jmckenzi/Downloads/PDM FILES/';
new = '/Users/jmckenzi/Desktop/PDM/';

% Find names of all raw files
rf = fileFinderAll(raw,'raw');
rf = rf(2:end,:);
numF = size(rf,1);

% Find names of the pdm files
pf = fileFinderAll(pdm,'pdm');

status = cell(numF,2);

% Loop through each raw file
for n = 1:numF
    
    % Update status
    status{n,1} = rf{n,2};
    
    % Need to ditch the 'Analyte...' part from the file... and add 'pdm'
    fx = strfind(rf{n,2},' Analyte')
    tmp = [rf{n,2}(1:fx-1) '.pdm'];
    
    fx = find(strcmp(pf(:,2),tmp));
    
    if numel(fx) == 1
        
        tmpn = [pf{fx,1} filesep pf{fx,2}];
        tmpl = [new pf{fx,2}];
        copyfile(tmpn,tmpl);
        
        status{n,2} = 'Pass';

    else
        status{n,2} = 'Fail';
    end
    
    
end



end


function exportFigure(name,type)
%exportFigure - save the current figure to the file name and format
%specified...

% There are various types possible...
switch type    
    case 'png'
        opt = '-dpng';
        
    case {'jpg','jpeg'}
        opt = '-djpeg';
        
    case 'eps'
        opt = '-depsc';
        
    otherwise
        % The default style...
        opt = 'dtiff';
        
end

% Fix the resolution
resol = '-r400';

% Alter page size
set(gcf, 'PaperPositionMode', 'auto');

% This function prints the figure
print(gcf,opt,resol,name);

    


end


function [ output_args ] = pubFig(name,type,resol)
%pubFig - save a publication ready figure, also in fig format
%
% By default save to Google Drive unless specified

% Depends on the computer
if ismac
    defP = '/Users/jmckenzi/Google Drive/Figures/';
else
    defP = 'E:\Drive\Figures\';
end

if ~exist(defP,'dir')
    defP = pwd;
end

% Check to see if name is a directory?
slash = strfind(name,filesep);
if numel(slash) > 0
    path = name(1:slash(end));
    
    if ~exist(path,'dir')
        % Use the default
        path = defP;
        name = [path name];
    end
else
    name = [defP name];
end

% Decide between type... PNG by default
switch type
    case 'eps'
        opt = '-depsc';        
    case {'jpeg','jpg'}
        opt = '-djpeg';        
    case {'tiff','tif'}
        opt = '-dtiff';        
    otherwise        
        opt = '-dpng';
end

% Resolution
if nargin == 3
    if isnumeric(resol)
        resol = ['-r' int2str(resol)];
        flag = true;
    else
        flag = false;
        resol = '-r100';
    end
else
    flag = false;
    resol = '-r100';
end

% Set the size to match what is on the screen!
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, opt, resol, name);

if flag
    print(gcf, '-dpng', '-r100', [name '-LowRes']);
end

% Save the Matlab figure too
%savefig(name);

end


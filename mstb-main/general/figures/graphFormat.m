function [ output_args ] = graphFormat( xLab, yLab, xItal, yItal, saveName, saveType )
%graphFormat Call this function to format the graph to a consistent
%quality.  Colours are variable, but let's try to keep to blue, red,
%magenta, green etc...

% Just a straight save, without additional formatting
if nargin == 2
    % What is the resolution?    
    resol = '-r100';
    
    % Select the format that you want...
    switch yLab
        case 'eps'
            opt = '-depsc';
        case 'png'
            opt = '-dpng';
            resol = '-r400';
        case 'jpeg'
            opt = '-djpeg';
        otherwise 
            opt = '-dtiff';
            resol = '-r400';
    end
    
    try
        % Set the size to match what is on the screen!
        set(gcf, 'PaperPositionMode', 'auto');
        print(gcf, opt, resol, xLab);
    
    catch error
        % Set the size to match what it on the screen!
        set(gcf, 'PaperPositionMode', 'auto');
        print(gcf, opt, resol, ['/Users/jmckenzi/Dropbox/Defaulted-' datestr(now, 'mm-dd') '-' datestr(now, 'HHMM')]);
        disp('There was an issue saving your file, but it is in dropbox somewhere');
        error
    end
    return
end

    

% Change the font?
%,'FontName', 'Garamond'

% Set the X and Y axes labels
switch xItal
    case 1
        xlabel(xLab, 'FontSize', 14, 'FontWeight', 'b', 'FontAngle', 'i');
    otherwise
        xlabel(xLab, 'FontSize', 14, 'FontWeight', 'b');
end

switch yItal
    case 1        
        ylabel(yLab, 'FontSize', 14, 'FontWeight', 'b', 'FontAngle', 'i');
    otherwise
        ylabel(yLab, 'FontSize', 14, 'FontWeight', 'b');
end

% Turn on the box!
box on;

% Legend?
lgnd = findobj(gcf, 'Tag', 'legend');
%legend(lgnd, 'Box', 'Off');
%leg = get(gcf, 'Children')

% Set the axes font sizes
set(gca, 'FontSize', 10, 'FontWeight', 'bold');

% Set the size to match what it on the screen!
%set(gcf, 'PaperPositionMode', 'auto');

if strcmp(saveName, '')
    % do nothing!
else
    
    % What is the resolution?
    resol = '-r300';
    
    % Select the format that you want...
    switch saveType
        case 'eps'
            opt = '-depsc';
        case 'png'
            opt = '-dpng';
            resol = '-r400';
        case 'jpeg'
            opt = '-djpeg';
        otherwise 
            opt = '-dtiff';
            resol = '-r400';
    end
    
    try
        % Set the size to match what it on the screen!
        set(gcf, 'PaperPositionMode', 'auto');
        print(gcf, opt, resol, saveName);
    
    catch error
        % Set the size to match what it on the screen!
        set(gcf, 'PaperPositionMode', 'auto');
        print(gcf, opt, resol, ['/Users/JSM/Dropbox/Defaulted-' datestr(now, 'mm-dd') '-' datestr(now, 'HHMM')]);
        disp('There was an issue saving your file, but it is in dropbox somewhere');
        error
end

end


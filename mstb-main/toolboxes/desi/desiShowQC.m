function desiShowQC(src,~,fig)
%desiShowQC - draw the plot for the file, if it exists

% Gett he guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Where to get the icons from?
dirIcon = deSlash([pwd filesep '/toolboxes/icons/']);

% Check to see if it all exists...
if ~isfield(dpn.opts,'qc')
    
    set(src,'CData',getImage(dirIcon,'gavel-2-48-warning'));
    hh = warndlg('No QC data available for this file');
    pause(1);
    delete(hh);
    return
end

% Now just run the function to draw the figure
[pass] = desiQCimages(dpn.opts.qc,true);

% Simple message...
if pass
    
    % Change the image to this...
    img = 'gavel-2-48-pass';
    
    hh = warndlg('PASS!','QC');
    pause(1);
    delete(hh);
else    
    
    % Set the image to this
    img = 'gavel-2-48-fail';
    
    errordlg('This file fails the QC test','QC','modal');
end

% Update the icon to show pass or fail
set(src,'CData',getImage(dirIcon,img));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ico] = getImage(path,name)
% Get the necessary icon.  Need to beautify it a little bit more...

[ico] = importdata([path name '.png']);
[ico] = iconProcess(ico.alpha,ico.cdata);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

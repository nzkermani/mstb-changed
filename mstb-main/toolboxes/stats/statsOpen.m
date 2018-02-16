function statsOpen(~,~,fig,defP)
%statsOpen - load a previously saved file

% THis is better than using a dialog box
wb = waitbar(1/2,'Opening a file');


% Ask for a file...
[a,b,c] = uigetfile({'*.mat'},'Load a File',defP);

% Continue?
if c == 0
    delete(wb);
    return
end

% Read in the sts variable of this file. What if it doesn't exist?
tmp = load([b a],'sts');
    
if ~isfield(tmp,'sts')
    warndlg('Invalid file - cannot find right bits');
    delete(wb);
    return
end

% Just clear all axes
statsAxesClearAll([],[],fig);

% Add the figure data
tmp.sts.fig = fig;

% Check that the magical 'obsID' exists. If not create it
if ~isfield(tmp.sts.raw,'obsID')
    vec = 1:size(tmp.sts.raw.sp,1);
    tmp.sts.raw.obsID = vec';
    
    % What about proc?
    procSz = size(tmp.sts.proc.sp,1);
    if procSz == size(tmp.sts.raw.sp,1)
        tmp.sts.proc.obsID = tmp.sts.raw.obsID;
    else   
        tmp.sts.proc.obsID = NaN(size(tmp.sts.proc.sp,1),1);
    end
end

% Save as guidata
guidata(fig.fig,tmp.sts);

% Now create a function that plots the data in the spectral window
statsPlotSpectra([],[],fig,'proc');

% Populate the table...
statsTablePopulate([],[],fig);

% Delete the waitbar
delete(wb);



end


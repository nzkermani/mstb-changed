function dpnSave(~,~,fig)
%dpnSave - save the current progress

% Guidata gather
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Strip out the graphics handles as they are not needed!
dpn.fig = [];

% THis is better than using a dialog box
wb = waitbar(1/2,'Saving');

% This part is taken from the matlab save function and just allows us to
% use the -v7.3 flag as default for large file sizes...
filters = {'*.mat','MAT-files (*.mat)'};

% Determine default file names etc...
try
    defP = dpn.file.dir;
    dots = strfind(dpn.file.nam,'.');
    defN = dpn.file.nam(1:dots(end)-1);
catch
    defP = dpn.imz.dir;
    dots = strfind(dpn.imz.nam,'.');
    defN = dpn.imz.nam(1:dots(end)-1);
end


% Just ask the user for a file name and location...
[fn,pn,~] = uiputfile(filters,...
    getString(message('MATLAB:uistring:filedialogs:SaveWorkspaceVariables')),...
    [defP defN '.mat']);

% If cancelled then nothing happens...
if ~isequal(fn,0)
    
    % Create the full length file name
    fn = strrep(fullfile(pn,fn), '''', '''''');
    
    % Save the file name - need to test that 'normal' files are ok with
    % this too.
    save(fn,'dpn','-v7.3');
    
end
    

% Now ask the user for a save location...
%uisave('dpn',[dpn.file.dir dpn.file.nam(1:end-6) '.mat']);

% Delete the waitbar
delete(wb);

end


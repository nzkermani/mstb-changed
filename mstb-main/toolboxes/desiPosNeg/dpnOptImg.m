function dpnOptImg(~,~,fig)
%dpnOptImg - upload an optical image or use one of the ion images...

% Gather the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return;
end

% Ask the user...
switch dpn.mode
    case 'dual'
        choice = questdlg('Chose an optical image...',...
            'Optical Image',...
            'MS1','MS2','Opt','Opt');
        
    case 'single'
        choice = questdlg('Chose an optical image...',...
            'Optical Image',...
            'MS1','Opt','Opt');
        
    otherwise
        error('Do not know what happened');
end

% Upload either optical or the MS image
switch choice
    
    % Upload the optical image as MS1
    case 'MS1'
        try
            img = repmat(dpn.d1.img,[1 1 3]);
            img = img - min(img(:));
            img = img ./ max(img(:));
            img = 1 - img;
        catch
            img = nansum(dpn.d1.sp,3);
        end
        
    % Upload the optical image to be MS2
    case 'MS2'
        try
            img = repmat(dpn.d2.img,[1 1 3]);
            img = img - min(img(:));
            img = img ./ max(img(:));
            img = 1 - img;
        catch
            img = nansum(dpn.d2.sp,3);
        end
        
    % Ask the user to provide an optical image
    case 'Opt'
        [imgName,imgPath] = uigetfile({'*.tif; *.png; *.tiff; *.jpg'},...
            'Select an optical image',...
            dpn.file.dir);
        
        if length(imgName) == 1
            return;
        end
        
        % Now read in the image...
        img = imread([imgPath imgName]);
       
        
    % There is no optical image
    otherwise        
        disp('abort');
        return
end

% Need to clear the existing axes to remove annotations etc...
if isfield(dpn,'anno')    
    f0 = findobj('Type','patch');
    delete(f0);        
    dpn = rmfield(dpn,'anno');
end

% Also consider removing existing coregistrations and the grid...
if isfield(dpn,'opt')
    dpn = rmfield(dpn,'opt');
    
    try
        newfig = rmfield(dpn.fig,'grid');
        dpn.fig = newfig;
    catch
        disp('no grid to delete');
    end
end

% Ensure that coregistration is disabled if we have picked to use the MS
% image instead of an optical image
switch choice(1:2)
    case 'MS'
        set(dpn.fig.tb.coreg,'Enable','off');
        set(dpn.fig.tb.fiducial,'Enable','off');
        dpn.opt.coreg = img;
    case 'Op'
        set(dpn.fig.tb.coreg,'Enable','on');
        set(dpn.fig.tb.fiducial,'Enable','on');
end


% Save the optical image to the structure
dpn.opt.orig = img;

% Update the guidata
guidata(fig.fig,dpn);

% Now that we have an image then plot it...
set(fig.ax.opt(2),'CData',img);
set(fig.ax.opt(1),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(img,2)],...
    'YLim',[1 size(img,1)]);

end


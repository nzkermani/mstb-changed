function desiRGB2Opt(src,event,fig)
%desiRGB2Opt - convert the RGB image to the optical image rather than
%uploading an H&E image

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Gather the image from the axes - need to ensure that there actually is
% some data here...
f0 = get(fig.ax.mv,'Children')
f1 = get(f0,'Type')
chk = strcmp(f1,'image')

fx = f0(chk)

if numel(fx) == 0
    disp('No image to be found. Bad luck');
elseif numel(fx) > 1
    disp('Error of unkonwn proportion');
else
    % Then there is just the one image, so use this!
    img = get(fx,'CData');
    
    % Need to clear the existing axes to remove annotations etc...
    if isfield(dpn,'anno')    
        f0 = findobj('Type','patch');
        delete(f0);        
        dpn = rmfield(dpn,'anno');
    end

    % Ensure that coregistration is disabled if we have picked to use the
    % MS image instead of an optical image
    set(dpn.fig.tb.coreg,'Enable','off');
    dpn.opt.coreg = img;
    
    % Save the optical image to the structure
    dpn.opt.orig = img;
    guidata(fig.fig,dpn);

    % Now that we have an image then plot it...
    set(fig.ax.opt(2),'CData',img);
    set(fig.ax.opt(1),'XTick',[],...
        'YTick',[],...
        'XLim',[1 size(img,2)],...
        'YLim',[1 size(img,1)]);

end




end


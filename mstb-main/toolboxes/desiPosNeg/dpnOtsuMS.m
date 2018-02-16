function dpnOtsuMS(src,event,fig,cor)
%dpnOtsuMS - do Otsu thresholding on the MS image

% Get the guidata
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Determine which mode is to be thresholded
tmp = get(cor.listMS,'String');
mode = get(cor.listMS,'Value');
mode = tmp{mode};
switch mode
    case 'Positive'
        dx = dpn.d1;
        ax = [dpn.fig.ax.ms1 dpn.fig.ax.ms2];
    case 'Negative'
        dx = dpn.d2;
        ax = [dpn.fig.ax.ms2 dpn.fig.ax.ms1];
    
    % These two options are those when just doing a single ion mode
    case 'Original'
        set(dpn.fig.ax.ms1(2),'CData',dpn.d1.img);
        set(dpn.fig.ax.ms1(1),'XTick',[],...
            'YTick',[],...
            'XLim',[1 size(dpn.d1.img,2)],...
            'YLim',[1 size(dpn.d1.img,1)]);
        return

    case 'Otsu'
        dx = dpn.d1;
        ax = [dpn.fig.ax.ms1 dpn.fig.ax.ms1];
        
end

% What is the state of the slider at the moment?
val = round(get(cor.slideMS,'Value'));

% Run the Otsu tresholding method
%[dx.tobg] = pixTOBG(dx.img,val,[]);
[dx.tobg,~] = dpnTOBG(dx.img,[],val);

% Perhaps we can consider smoothing the image here...
filt = fspecial('average',3);
dx.tobg = filter2(filt,dx.tobg) >  0.75;

% Plot the tobg image
set(ax(4),'CData',dx.tobg);
set(ax(3),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(dx.tobg,2)],...
    'YLim',[1 size(dx.tobg,1)]);

% Update the main image
if isfield(dpn,'d2')
    set(ax(2),'CData',dx.img);
    set(ax(1),'XTick',[],...
        'YTick',[],...
        'XLim',[1 size(dx.img,2)],...
        'YLim',[1 size(dx.img,1)]);
end

% Save to the guidata
switch mode
    case 'Positive'
        dpn.d1 = dx;
    case 'Negative'
        dpn.d2 = dx;
    otherwise
        % For single ion mode
        dpn.d1 = dx;
end
guidata(fig.fig,dpn);


end


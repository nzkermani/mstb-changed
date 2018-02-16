function pdDisplay(~,~,fig)
%pdDisplay - show the images as requested

% Get the guidata
data = guidata(fig.fig);
if isempty(data)
    return
end

% Now work out the various options for visualisation
opts.ionVals = fig.ionDisplay.Value;
opts.doLog = fig.doLog.Value;
opts.imin = str2double(fig.minInt.String);
opts.imax = str2double(fig.maxInt.String);

% Now we need to generate the image...
img = nansum(data.image(:,:,opts.ionVals),3);

% Scale between 0 and 1
img = img / max(img(:));

% Set values outside of these values
img(img < opts.imin) = opts.imin;
img(img > opts.imax) = opts.imax;

% Log or not
if opts.doLog
    img = img * 1000;
    os = nanmedian(img(img > 0));
    img = log(img + os);
end

% Now we need to change the figure to match
set(fig.ax(2),'CData',img);

xlim(fig.ax(1),[0.5 size(img,2)+0.5]);
ylim(fig.ax(1),[0.5 size(img,1)+0.5]);

end


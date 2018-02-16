function [ flag ] = desiQCimages(qc,plot)
%desiQCimages - make a plot showing the QC images

% These are the raffinose peaks in neg and pos
raff = [503.1627 527.1583];

% Rather than plot both, we can assume that the most intense of the two
% ions is the correct one, hence telling us the ion mode.
sz = size(qc.sp);
avgs = nanmean(reshape(qc.sp,[sz(1) * sz(2) sz(3)]),1);
[~,i] = max(avgs);

% Determine ppm diffs...
ppm = 1e6 * (qc.mz(:,:,i) - raff(i)) ./ raff(i);


mz = qc.mz(:,:,i);
sp = qc.sp(:,:,i);
mask = abs(ppm) > 1;

flag = sum(mask(:)) / numel(mask);
if flag < 0.1
    flag = true;
else
    flag = false;
end

% No figure if not requested
if nargin == 1
    return    
end

% Draw a figure...
figure('Units','normalized',...
    'Position',[0 0 1 1]);

ax = zeros(1,6);

%%%%%%%%%% Exact
% NEG - image
ax(1) = subplot(2,3,1);
tmp = sp;
tmp(mask) = NaN;
imagesc(tmp);
title(['Raffinose, m/z = ' sprintf('%0.4f',raff(i)) ' +/- 1 ppm'],...
    'FontSize',14,'FontWeight','bold');
y = colorbar;
ylabel(y,'Raw intensity','FontSize',12,'FontWeight','bold');

% NEG - mz
ax(2) = subplot(2,3,2);
tmp = mz;
tmp(mask) = NaN;
imagesc(tmp);
title(['Raffinose, m/z = ' sprintf('%0.4f',raff(i)) ' +/- 1 ppm'],...
    'FontSize',14,'FontWeight','bold');
y = colorbar;
ylabel(y,'m/z','FontSize',12,'FontWeight','bold');

% NEG - ppm
ax(3) = subplot(2,3,3);
tmp = ppm;
tmp(mask) = NaN;
imagesc(tmp)
title(['Raffinose, m/z = ' sprintf('%0.4f',raff(i)) ' +/- 1 ppm'],...
    'FontSize',14,'FontWeight','bold');
y = colorbar;
ylabel(y,'\delta ppm','FontSize',12,'FontWeight','bold');

%%%%%%%%%% Total
% NEG - image
ax(4) = subplot(2,3,4);
imagesc(qc.sp(:,:,i));
title(['Closest Peak to m/z = ' sprintf('%0.4f',raff(i))],...
    'FontSize',14,'FontWeight','bold');
y = colorbar;
ylabel(y,'Raw intensity','FontSize',12,'FontWeight','bold');

% NEG - mz
ax(5) = subplot(2,3,5);
imagesc(qc.mz(:,:,i));
title(['Closest Peak to m/z = ' sprintf('%0.4f',raff(i))],...
    'FontSize',14,'FontWeight','bold');
y = colorbar;
ylabel(y,'m/z','FontSize',12,'FontWeight','bold');

% NEG - ppm
ax(6) = subplot(2,3,6);
imagesc(ppm)
title(['Closest Peak to m/z = ' sprintf('%0.4f',raff(i))],...
    'FontSize',14,'FontWeight','bold');
y = colorbar;
ylabel(y,'\delta ppm','FontSize',12,'FontWeight','bold');

set(ax,'XTick',[],'YTick',[]);

end


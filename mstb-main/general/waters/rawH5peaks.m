function [pks,sm] = rawH5peaks(mz,sp)
%rawH5peaks - perform peak picking on an average spectrum


% What about Savitzky-Golay filtering?
sm = sgolayfilt(sp,7,21);
%sm = sp;

% Can we define a prominence threshold with which to ditch low intensity
% peaks from the get go?
%minIntensity = max(sm) * 0.0005;
%minIntensity = nanmedian(sm) * 3;% * 5
minIntensity = median(sm) * 2;% * 5

% figure; hold on;
% plot(mz,sp,'red');
% plot(mz,sm,'blue');

figure;

findpeaks(sm,mz,...
    'MinPeakProminence',minIntensity,...
    'MinPeakWidth',0.01,...
    'MaxPeakWidth',0.1,...
    'Annotate','extents',...
    'WidthReference','halfprom');

[pks.ints,pks.mz,pks.span,pks.prom] = findpeaksJSM(sm,mz,...
    'MinPeakProminence',minIntensity,...
    'MinPeakWidth',0.01,...
    'MaxPeakWidth',0.1,...
    'Annotate','extents',...
    'WidthReference','halfprom');

end


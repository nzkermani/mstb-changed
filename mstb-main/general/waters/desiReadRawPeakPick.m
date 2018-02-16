function [ output_args ] = desiReadRawPeakPick(mz,sp,fq)
%desiReadRawPeakPick - determine peaks and stuff from waters mz/sp/fq
%information

% Need to perform a fair bit of smoothing of the sp in order to be able to
% determine local maxima
tic
[mean10] = movingWindow(sp,10,'mean');
toc
tic
[mean5] = movingWindow(sp,4,'mean');
toc

% Resmooth
rem10 = movingWindow(mean10,10,'mean');
rem5  = movingWindow(mean5,5,'mean');


figure; hold on;
plot(mz,sp/max(sp),'k');
plot(mz,mean10/max(mean10),'r');
plot(mz,mean5/max(mean5),'b');

plot(mz,rem10/max(rem10),'m');
plot(mz,rem5/max(rem5),'c');




end


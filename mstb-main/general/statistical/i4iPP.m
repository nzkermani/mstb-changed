function i4iPP(ii,ss,ss2)
%i4iPP - plot various graphs showing power and sample sizes etc...

% Find variables that have a q-value less than 0.01
fx = ss.pq(:,1) < 0.05;

% figure; hist(ss.pow(fx,1),0.75:0.01:1);
% set(gca,'FontSize',14);
% xlabel('Power','FontWeight','bold');
% ylabel('Frequency','FontWeight','bold');
% %xlim([0.78 1.01])

figure; hist(ss2.pow(fx,2),0:5:70);%,0:5:50);%0:5:140);
set(gca,'FontSize',14);
xlabel('Sample size','FontWeight','bold');
ylabel('Frequency','FontWeight','bold');
xlim([-5 70])

mean(ss2.pow(fx,2))
std(ss2.pow(fx,2))

return

avg = nanmean(ii.sp,1);


plotXYC(ii.mz',avg,-log10(ss.pq(:,1)'),'-log_1_0 p-value');

pp = ss.pow(:,2);
mask = pp > 500;
pp(mask) = 500;
plotXYC(ii.mz',avg,pp','Sample size');

mm = median(ss.pow(:,2))
sum(ss.pow(:,2) < mm)
median(ss.pow(~mask,2))
figure; hist(ss.pow(~mask,2),50);
set(gca,'FontSize',14);
xlabel('Sample size','FontWeight','bold');
ylabel('Frequency','FontWeight','bold');
figure; scatter(ss.pow(~mask,2),-log10(ss.pq(~mask,2)));

return

% Find me variables that are significant with q < 0.01
alpha = 1.9;
fx = ss.pq(:,2) < alpha;

% Find me variables with a sample size < 1000
fy = ss.pow(:,2) < 1000;

% COmbine various metrics together
fa = fx & fy;

figure; hist(ss.pow(fa,2),0:10:1000);


end


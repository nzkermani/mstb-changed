function [ output_args ] = dpnPaperSpectrumPlot(dpn,m1,m2,cols)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Define offsets
%os = [5 4 3 2 1 0];
os = [(size(m1,1)*2)-1:-1:0]

%mz = [600 1000];

figure('Units','pixels','Position',[781 652 949 424]); 
hold on;

% DS 1
os1 = os(1:size(m1,1));%os(1:2:end)
for n = 1:size(m1,1)
    
    sp = (m1(n,:)/max(m1(n,:))) * 0.95;
    [a,b] = insertZeros(dpn.d1.mz,sp,0.001);
    plot(a,b + os1(n),'Color',cols(n,:),...
        'LineWidth',2);

end

% DS 2
os2 = os(size(m2,1)+1:end);%os(2:2:end)
for n = 1:size(m2,1)
    
    sp = (m2(n,:)/max(m2(n,:))) * 0.95;
    [a,b] = insertZeros(dpn.d2.mz,sp,0.001);
    plot(a,b + os2(n),'Color',cols(n,:),...
        'LineWidth',2);

end

xlim([150 1000]);
ylim([-0.01 max(os) + 1.01]);

text(120,mean(os1)+0.5,'Positive Mode',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'FontSize',18);

text(120,mean(os2)+0.5,'Negative Mode',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'FontSize',18);

box on;

set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal',...
    'FontSize',16);

xlabel('m/z','FontSize',16,'FontWeight','bold');

grid off;


return

i = 1;
sp = (m1(1,:)/max(m1(1,:))) * 0.95;
[a,b] = insertZeros(dpn.d1.mz,sp,0.001);
plot(a,b + os(i),'Color',cols(1,:),...
    'LineWidth',2);

i = 2;
sp = (m1(2,:)/max(m1(2,:))) * 0.95;
[a,b] = insertZeros(dpn.d1.mz,sp,0.001);
plot(a,b + os(i),'Color',cols(2,:),...
    'LineWidth',2);

i = 3;
sp = (m2(1,:)/max(m2(1,:))) * 0.95;
[a,b] = insertZeros(dpn.d2.mz,sp,0.001);
a(end+1) = 1500;
b(end+1) = 0;
plot(a,b + os(i),'Color',cols(1,:),...
    'LineWidth',2);

i = 4;
sp = (m2(2,:)/max(m2(2,:))) * 0.95;
[a,b] = insertZeros(dpn.d2.mz,sp,0.001);
a(end+1) = 1500;
b(end+1) = 0;
plot(a,b + os(i),'Color',cols(2,:),...
    'LineWidth',2);

xlim([150 1000]);
ylim([-0.01 max(os) + 1.01]);

text(120,3,'Positive Mode',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'FontSize',18);

text(120,1,'Negative Mode',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'FontSize',18);

box on;

set(gca,'YTickLabel',[],...
    'YTick',0,...
    'XTickMode','auto',...
    'XTickLabelMode','auto',...
    'LineWidth',5,...
    'TickLength',[0 0],...
    'YDir','normal',...
    'FontSize',16);

xlabel('m/z','FontSize',16,'FontWeight','bold');

grid off;


end


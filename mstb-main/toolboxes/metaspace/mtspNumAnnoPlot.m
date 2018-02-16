function [ output_args ] = mtspNumAnnoPlot(fName,numA)
%mtspNumAnnoPlot - plot the annotation quantities

chk = max(numA,[],2) == 0;

fName = fName(~chk,:);
numA = numA(~chk,:);


fx = numA(:,1) > numA(:,2);


figure; hold on;

scatter(numA(~fx,1),numA(~fx,2),80,'blue','o','filled');

scatter(numA(fx,1),numA(fx,2),100,'red','o','filled');

llim = max(numA(:)) * 1.1;

line([0 llim],[0 llim],'Color','k',...
    'LineWidth',2);

xlim([0 llim]);
ylim([0 llim]);


% Add text elements for the bad files
idx = find(fx);
for n = 1:numel(idx)
    
    text(numA(idx(n),1)+5,numA(idx(n),2),fName{idx(n),1},...
        'FontSize',14);
    
end

box on;
set(gca,'FontSize',14);

xlabel('Annotations (original)','FontWeight','bold','FontSize',16);
ylabel('Annotations (recalibrated)','FontWeight','bold','FontSize',16);



end


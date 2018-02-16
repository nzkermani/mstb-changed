function statsDrawROC(roc)
%statsDrawROC - make a pretty curve

%[oNO] = findOptimal(roc.x,roc.yNoCV(2:end));
%[oCV] = findOptimal(roc.x,roc.yCV(2:end));

%return

% Determine the x values that we should be using...
if numel(roc.x) ~= numel(roc.yNoCV)
    x1 = roc.xNoCV;
    y1 = roc.yNoCV;
else
    x1 = roc.x;
    y1 = roc.yNoCV(2:end);
end

if numel(roc.x) ~= numel(roc.yCV)
    x2 = roc.xCV;
    y2 = roc.yCV;
else
    x2 = roc.x;
    y2 = roc.yCV(2:end);
end

figure('Units','pixels','Position',[965 750 712 588]); hold on;
ax = zeros(1,2);

ax(1) = plot(x1,y1,'r',...
    'LineWidth',2);

ax(2) = plot(x2,y2,'b',...
    'LineWidth',2);

line([0 1],[0 1],'Color','k','LineStyle','--');
%line([0 0.5],[1 0.5],'Color','k','LineStyle','--');


scatter(0.3,0.2,200,'r','o','filled');
scatter(0.3,0.15,200,'b','o','filled');

text(0.35,0.2,['AUC = ' sprintf('%0.3f',roc.aucNoCV) ' (not cross validated)'],...
    'FontSize',18);

text(0.35,0.15,['AUC = ' sprintf('%0.3f',roc.aucCV) ' (cross validated)'],...
    'FontSize',18);


% legend(ax,{'ROC - Not Cross Validated';'ROC - Cross Validated'},...
%     'Location','southeast');


xlabel('FDR','FontWeight','bold');
ylabel('TPR','FontWeight','bold');
%scatter(oNO(1),oNO(2),80,'r','o','filled','MarkerEdgeColor','k');
%scatter(oCV(1),oCV(2),80,'b','o','filled','MarkerEdgeColor','k');

box on;

axis square

set(gca,'FontSize',18);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [o] = findOptimal(x,y)
% Determine the point closest to [0,1]

x = x;
y = 1 - y;

figure; stem(x,x+y);




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function annotateClassInfo(lm,ass)
%annotateClassInfo - tell us what kind of annotations are predominant.
%
% INPUTs
% lm    - database
% ass   - ass structure from annotateMZ.m

% Mask...
mask = ass.dbMatch;

makeGraph(lm.Class2,mask);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = makeGraph(cl,mask)
% Given a single class make a graph showing composition

% What about all lipids in the database?
[all,~,allind] = unique(cl);
numA = numel(all);
freq = zeros(numA,4);

% Loop through
for n = 1:numA
    
    % Indices of this class
    fx = allind == n;
    
    % Indices of those that were annotated
    fy = fx(mask);
    
    % Save...
    freq(n,1:2) = [sum(fy) sum(fx)];
    
end

% Save proportions
freq(:,3) = freq(:,1) / sum(freq(:,1));
freq(:,4) = freq(:,1) ./ freq(:,2);

% Strip out all zeros
fx = freq(:,1) == 0;
freq(fx,:) = [];
all(fx,:) = [];
numA = numA - sum(fx);

% Sort the data...
[~,idx] = sortrows(freq,-1);
freq = freq(idx,:);
all = all(idx);

% Determine colours
try
    cols = parula(numA);
catch
    cols = jet(numA);
end

% Draw a figure
figure('Units','normalized','Position',[0 0.5 1 0.4],'Color',[0.95 0.95 0.95]); 

subplot(1,3,1);
hold on;
for n = 1:numA
    bar(n,freq(n,1),'FaceColor',cols(n,:));
    text(n,freq(n,1)+2,all{n},'Rotation',30,'FontSize',14);
end
set(gca,'XTick',[],'XTickLabels',[],'Color',[0.95 0.95 0.95],'FontSize',14);
xlim([0.5 numA+0.5]);
ylim([0 max(freq(:,1)) * 1.2]);
ylabel('Number of annotations','FontSize',18,'FontWeight','bold');

% Pie chart - simple one
subplot(1,3,2);
ff = pie(freq(:,3));%,expl);
try
    set(ff(2:2:end),'FontSize',14);
catch
    set(ff(1:2:end),'FontSize',14);
end

% Final bar chart showing proportion of annotated classes
subplot(1,3,3);
hold on;
for n = 1:numA
    bar(n,freq(n,4)*100,'FaceColor',cols(n,:));
end
set(gca,'XTick',[],'XTickLabels',[],'Color',[0.95 0.95 0.95],'FontSize',14);
xlim([0.5 numA+0.5]);
%ylim([0 100]);
ylabel('Percentage of annotations','FontSize',18,'FontWeight','bold');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
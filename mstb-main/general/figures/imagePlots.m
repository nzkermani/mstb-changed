function [ output_args ] = imagePlots( tr,te )
%imagePlots - all sorts of plots for Jocelyn

gridPlot(tr,te,'prob')

scoresScatter(tr,te)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gridPlot(tr,te,type)


fig.fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);
%hold on;

ax0 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 0.5 1]);

posns = [0.5 0.5 0.25 0.5;...
    0.75 0.5 0.25 0.5;...
    0.5 0 0.25 0.5;...
    0.75 0 0.25 0.5];
ax = zeros(4,1);
for n = 1:4
    ax(n,1) = axes('Parent',fig.fig,...
        'Units','normalized',...
        'Position',posns(n,:));
end


% Plots...
switch type
    case 'prob'
        tmp = tr.prob(:,:,1:3);        
    case 'score'
        tmp = imScale(tr.score(:,:,1:3));
                
    otherwise
        if isnumeric(type)
            tmp = tr.prob(:,:,1:3) > type;
            type = 'thresh';
        else
            type = prob;
            tmp = tr.prob(:,:,1:3);
        end
end

axes(ax0);
imagesc(tmp);
text(1,5,tr.name,...
    'FontSize',16,...
    'FontWeight','bold',...
    'Color','white');

% Now for the subplots...
for n = 1:4
    
    axes(ax(n,1));
    
    switch type
        case 'prob'
            imagesc(te(n).prob(:,:,1:3));
            
        case 'score'
            imagesc(imScale(te(n).score(:,:,1:3)));
            
        case 'thresh'
            imagesc(te(n).score(:,:,1:3) > thresh);
            
        otherwise
            imagesc(te(n).score(:,:,1:3));
    end
    
    text(1,5,te(n).name,...
        'FontSize',16,...
        'FontWeight','bold',...
        'Color','white');
end
        
            





end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scoresScatter(tr,te)

figure; hold on;

% Scatter plot of the train set scores
t1 = tr.score(:,:,1); 
t1 = t1(:);
t2 = tr.score(:,:,2); 
t2 = t2(:);
idx = randperm(numel(t1),600);

scatter(t1(idx),t2(idx),80,'k','o','filled');

cols = jet(size(te,2));
for n = 1:size(te,2)
        
    u1 = te(n).score(:,:,1);
    u1 = u1(:);
    u2 = te(n).score(:,:,2);
    u2 = u2(:);
    idx2 = randperm(numel(u1),600);

    scatter(u1(idx2),u2(idx2),50,cols(n,:),'d','filled',...
        'MarkerEdgeColor','k');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
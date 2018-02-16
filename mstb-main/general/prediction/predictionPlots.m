function predictionPlots(data)
%imagePlots - all sorts of plots for Jocelyn

gridPlot(data,'thresh');
%gridPlot(data,'score');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gridPlot(data,type)

thresh = 0.1;

fig.fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);

posns = [0 0 0.5 1;...
    0.5 0.5 0.25 0.5;...
    0.75 0.5 0.25 0.5;...
    0.5 0 0.25 0.5;...
    0.75 0 0.25 0.5];
numAx = size(posns,1);

ax = zeros(numAx,1);

for n = 1:numAx

    ax(n,1) = axes('Parent',fig.fig,...
        'Units','normalized',...
        'Position',posns(n,:));
end

% Now for the subplots...
for n = 1:numAx
    
    axes(ax(n,1));
    
    % Apply the tobg to make background pixels gray, i.e. excluded from
    % analysis...
    if size(data(n).tobg,2) == 1
        tobg = reshape(data(n).tobg,data(n).sz(1:2));
    end
    
    switch type
        case 'prob'
            tmp = data(n).prob(:,:,1:3);
            for r = 1:3
                ff = tmp(:,:,r);
                ff(~tobg) = 1;
                tmp(:,:,r) = ff;
            end
            imagesc(tmp);
            
        case 'score'
            tmp = data(n).prob(:,:,1:3);
            for r = 1:3
                ff = tmp(:,:,r);
                ff(~tobg) = 1;
                tmp(:,:,r) = ff;
            end
            imagesc(tmp);
            
        case 'thresh'
            tmp = data(n).prob(:,:,1:3);
            for r = 1:3
                ff = tmp(:,:,r);
                ff(~tobg) = 0;
                tmp(:,:,r) = ff;
            end
            imagesc(tmp > thresh);
            
        otherwise
            imagesc(data(n).score(:,:,1:3));
    end
    
    set(gca,'XTick',[],'YTick',[]);
    
    text(1,5,data(n).name,...
        'FontSize',16,...
        'FontWeight','bold',...
        'Color','white');
end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scoresScatter(data)

figure; hold on;

% Scatter plot of the scores
cols = jet(size(data,2));
for n = 1:size(data,2)
        
    u1 = data(n).score(:,:,1);
    u1 = u1(:);
    u2 = data(n).score(:,:,2);
    u2 = u2(:);
    idx2 = randperm(numel(u1),600);

    scatter(u1(idx2),u2(idx2),50,cols(n,:),'d','filled',...
        'MarkerEdgeColor','k');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
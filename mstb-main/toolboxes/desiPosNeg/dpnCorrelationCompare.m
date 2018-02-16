function dpnCorrelationCompare(op)
%dpnCorrelationCompare - plots to explain why there is one pixel which
%doesn't look very good

% Which figure to produce?
figNum = 1;
switch figNum
    
    case 1
        xy = [59 36];

        % These are the three problematic variables
        mzs = [311.1695 325.1842 339.2000]; %297.1537 
        mzText = {'311.170';'325.184';'339.200'};

    case 2
        xy = [32 40];

    case 3
        xy = [56 8];
        mzs = [255.2340 281.2496];
        mzText = {'255.234';'281.250'};

end

% Draw the figure
figure('Units','pixels','Position',[588 497 1633 691]); 

% Add the parts
ax = subplot(1,2,1); hold on;
scatter([0.5 xy(2)],[xy(1) 0.5],80,'k','d','filled');
set(gca,'YDir','reverse');

imagesc(op.crr);
axis tight
axis square
axis off
cmap = hot(100);
colormap(cmap);
cb = colorbar;
set(gca,'FontSize',16);
ylabel(cb,'Correlation Coefficient','FontSize',18,'FontWeight','bold');
set(cb,'FontSize',16);


%mzF = find(mzFind(op.mz1,mzs,4));

% X and Y intensities
x = squeeze(op.n1(xy(1),xy(2),:));
y = squeeze(op.n2(xy(1),xy(2),:));

% ...but which of these are interpolated intensities?
i1 = op.interp1(xy(1),xy(2));
i2 = op.interp2(xy(1),xy(2));

% CHange intensities for scatter plot if required...
if i2 == 1 && i1 == 0
    y2 = x;
    x = y;
    y = y2;
end


% Here we need to find the dodgy ones?
switch figNum    
    case 1
        mask = x < 500 & y > 1250;
        ofset = [100 -50; 100 0; 100 +50];
    case 2
        mask = x < 500 & y > 1250;
    case 3
        mask = x > 3000;
        ofset = [50 200; 50 200];
end
        
x(mask)
y(mask)
op.mz1(mask)

corr(x,y)
corr(x(~mask),y(~mask))

ax(2) = subplot(1,2,2); hold on;

scatter(x,y,...
    70,'k','o','filled',...
    'MarkerEdgeColor',[0.8 0.8 0.8]);

% Only add predefined bits and pieces
if sum(mask) > 0
    annos = [x(mask) y(mask)];
    locn = annos + ofset;
    %locn(:,1) = 300;

    text(locn(:,1),locn(:,2),mzText,...
        'FontSize',16);

    for n = 1:numel(mzText)
        line([locn(n,1) annos(n,1)],[locn(n,2) annos(n,2)])
    end

    scatter(x(mask),y(mask),...
        70,'r','o','filled',...
        'MarkerEdgeColor','k');
end

% What about a line?
xl = xlim
yl = ylim
mxx = max([xl yl])
line([0 mxx],[0 mxx],'LineStyle','--');


% Now just format as required
box on;
xlabel('Intensities: Interpolated','FontSize',16);
ylabel('Intensities: Measured','FontSize',16);
set(gca,'FontSize',16);

return

figure; hold on;
stem(squeeze(op.n1(:,xy(2)-1:xy(2)+1,mzF(1))));
stem(squeeze(-op.n2(:,xy(2)-1:xy(2)+1,mzF(1))));

figure;
a(1) = subplot(1,2,1); imagesc(op.n1(:,:,mzF(1)));
a(2) = subplot(1,2,2); imagesc(op.n2(:,:,mzF(1)));
linkaxes(a,'xy');
end


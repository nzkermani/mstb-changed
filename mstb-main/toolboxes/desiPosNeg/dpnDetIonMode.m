function [ion,fig,flag] = dpnDetIonMode(sp)
% Determine the ion modes of the various lines, considering that they
% haven't all been done perfectly.
%
% Tested using (JLA003_tumour_S4)
%
%
%
% James McKenzie, 2016.




% Note that mzcounts is just the number of features detected, rather than
% the intensities. If we calculate another metric perhaps it will be easier
spsum = zeros(size(sp.mzcounts));
spmax = zeros(size(sp.mzcounts));
for r = 1:size(sp.X,1)
    px = sp.pixel(r,1);
    py = sp.pixel(r,2);
    spsum(px,py) = sum(sp.counts{r,1});
    spmax(px,py) = max(sp.counts{r,1});
    
end
rowSum = sum(spsum,1)';
rowMax = max(spmax,[],1)';

% Maybe normalise to scale between 0 and 1?
rowSum = rowSum / max(rowSum);
rowMax = rowMax / max(rowMax);


% First of all we run the function with just 2 clusters.  If this goes
% wrong then just use the maximum rather than mean and max
[ion,fig,flag] = doDetermination(sp,rowSum,rowMax,2);

% if flag
%     [ion,fig,flag] = doDetermination(sp,rowMax,ones(size(rowSum)),2);
% end    

% Now if this is also flagged as being a poor determination of ion modes,
% then we need to draw a GUI for the user to specify which scans belong to
% which modes
if flag
    % Draw the GUI here...
    [f2] = drawDPNgui(sp);
    uiwait(f2.fig);
    try
        ion = 1 + cell2mat(get(f2.bg(:,2),'Value'));
    catch
        disp('User closed window before determination');
        flag = true;
        return
    end
       
    
    % These are the proper outputs as expected...
    %fig = fig.fig;
    flag = false;
    
    % Close the figure that we used to determine scan ownership
    close(f2.fig);
    
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ion,fig,flag] = doDetermination(sp,rowSum,rowMax,numClust)

% Simple calculations...
nRow = size(sp.mzcounts,2);
xvals = 1:nRow;

% Clustering?
y = pdist([rowSum rowMax],'euclidean');
z = linkage(y,'complete');
fx = cluster(z,'maxclust',numClust);
%fx = cluster([rowSum rowMax],'maxclust',numClust,'Distance','seuclidean');
fxHi = fx == 1;
fxLo = fx == 2;

% Determine if there are any 'three in a row'. If so, then consider that
% this is an error...
[flag] = threeRow(fx);

% Draw the figure
totSum = reshape(sp.X,size(sp.mzcounts));
%totSum = nansum(totSum,1);

[fig] = makeFigure2(totSum,fxHi,fxLo,nansum(totSum,1),xvals,flag);
%[fig] = makeFigure(sp,fxHi,fxLo,rowSum,xvals,flag);

% Extract the 1/2 matrix as expected
[ion] = oneortwo(fxHi);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = makeFigure(sp,fxHi,fxLo,rowSum,xvals,flag)

% Draw a new figure
fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);
hold on;

% Full image
ax(1) = subplot(2,2,1);
imagesc(sp.mzcounts);
title('Full image','FontSize',20);

% Spectral sums
ax(2) = subplot(2,2,3); hold on;
plot(xvals(fxHi),rowSum(fxHi),'-or','MarkerFaceColor','r','MarkerSize',10);
stem(xvals(fxHi),rowSum(fxHi),':or');
plot(xvals(fxLo),rowSum(fxLo),'-ob','MarkerFaceColor','b','MarkerSize',10);
stem(xvals(fxLo),rowSum(fxLo),':ob');
xlim([0.5 numel(xvals)+0.5]);
title('Scan intensities','FontSize',20);

% High intensity image...
ax(3) = subplot(2,2,2); hold on;
imagesc(log(sp.mzcounts(:,fxHi)));
scatter(1,1,200,'r','o','filled','MarkerEdgeColor','w');
xlim([0.5 sum(fxHi)+0.5]);
ylim([0.5 size(sp.mzcounts,1)+0.5]);
set(gca,'YDir','reverse');
title('Image part 1','FontSize',20);

% Low intensity image...
ax(4) = subplot(2,2,4); hold on;
imagesc(log(sp.mzcounts(:,fxLo)));
scatter(1,1,200,'b','o','filled','MarkerEdgeColor','w');
xlim([0.5 sum(fxHi)+0.5]);
ylim([0.5 size(sp.mzcounts,1)+0.5]);
set(gca,'YDir','reverse');
title('Image part 2','FontSize',20);

colormap(redbluecmap);

% Draw something to say that the splitting of the modes is probably not
% optimal
if flag
    annotation('textbox',[0.1 0.425 0.8 0.2],...
        'String','>>> Fail <<<',...
        'FontSize',40,...
        'Color','red',...
        'EdgeColor','none',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = makeFigure2(sp,fxHi,fxLo,rowSum,xvals,flag)

% Draw a new figure
fig = figure('Units','normalized',...
    'Position',[0.25 0.1 0.3 0.9]);

hAxis = axes('units','normalized','pos',[0 0 1 1],...
    'visible','off','handlevisibility','off');

text(0.025,0.96,'a)',...
    'Parent',hAxis,...
    'FontSize',20,...
    'FontWeight','bold');

text(0.025,0.5,'b)',...
    'Parent',hAxis,...
    'FontSize',20,...
    'FontWeight','bold');


ax = [0 0 0];

% Full image
ax(1) = axes('Parent',fig,...
    'Units','normalized',...
    'Position',[0.09 0.525 0.9 0.45]);
imagesc(log(sp));
colormap(gray);
set(ax(1),'XTick',[],'YTick',[]);
%title('Full image','FontSize',20);

% Spectral sums
ax(2) = axes('Parent',fig,...
    'Units','normalized',...
    'Position',[0.09 0.065 0.9 0.45]);
%stem(xvals,rowSum);
hold on;
plot(xvals(fxHi),rowSum(fxHi),'-or',...
    'MarkerFaceColor','r','MarkerSize',8,...
    'MarkerEdgeColor',[0.8 0.8 0.8]);
stem(xvals(fxHi),rowSum(fxHi),':or');
plot(xvals(fxLo),rowSum(fxLo),'-ob',...
    'MarkerFaceColor','b','MarkerSize',8,...
    'MarkerEdgeColor',[0.8 0.8 0.8]);
stem(xvals(fxLo),rowSum(fxLo),':ob');

box on;
xlim([0.5 numel(xvals)+0.5]);
%ylim([0 1.05]);
set(ax(2),'XTick',[],'YTick',[]);
xlabel('Scan','FontSize',18,'FontWeight','bold');
ylabel('Scan intensity (sum)','FontSize',18,'FontWeight','bold');

%title('Scan intensities','FontSize',20);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ion] = oneortwo(fx)
% Make a vector of 1s and 2s to identify the scans...

ion = double(fx) + 1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = threeRow(fx)
% Determine if there are three in a row

flag = false;
for n = 1:numel(fx)-2
    
    if fx(n) == fx(n+1) && fx(n) == fx(n+2)
        flag = true;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawDPNgui(sp)

fig.fig = figure('Units','normalized',...
    'Position',[0.1 0.25 0.8 0.5]);

fig.ax = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0.3 1 0.7]);

imagesc(sp.mzcounts,'Parent',fig.ax);
colormap(redbluecmap);
caxis(prctile(sp.mzcounts(:),[10 90]));

fig.pan = uipanel('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0.1 1 0.2]);

numS = size(sp.mzcounts,2);

wd = 1 / numS;

fig.bg = zeros(numS,3);

pos = 0;
for n = 1:numS
    
    fig.bg(n,1) = uibuttongroup('Parent',fig.pan,...
        'Units','normalized',...
        'Position',[pos 0 wd 1]);
    
    % Add the two modes
    for r = 1:2
        if r == 1
            ppos = [0 0 1 0.5];
        else
            ppos = [0 0.5 1 0.5];
        end
        
        if mod(n,2) == 0
            val = 1;
        else
            val = 0;
        end
        
        fig.bg(n,r+1) = uicontrol('Parent',fig.bg(n,1),...
            'Style','radiobutton',...
            'Units','normalized',...
            'Position',ppos,...
            'Value',val);
    end
            
    pos = pos + wd;
    
end
    
% Draw a button to finish the thing
fig.finish = uicontrol('Parent',fig.fig,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.25 0 0.5 0.05],...
    'String','FINISH!',...
    'Callback','uiresume(gcf)');

disp('Ready to go');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
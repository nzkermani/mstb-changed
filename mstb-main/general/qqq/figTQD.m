function figTQD(all)
%figTQD - make a fancy figure

% How many channels?
numC = size(all,2);

% Determine the maximum values
%maxVals = max(reshape(img,[size(img,1)*size(img,2) size(img,3)]),[],1);

f0 = findobj('Name','TQ-DESI'); close(f0);

% del2 = [12 23 34 45 56 67 78 89]; 
% for n = 1:size(all,2)
%     
%     all(n).img(:,del2,:) = [];
%         
% end

% Determine background pixels using the third image...
sz = all(1).sz;
bg = zeros(sz(1:2));
for n = 1:numC
    bg = bg + nansum(all(n).img,3);
end
bg = reshape(bg,[sz(1)*sz(2) 1]);
bgPix = getObjectPixels(bg);
bgPix = reshape(bgPix,[sz(1) sz(2)]);

% Smooth the image somewhat
fh = fspecial('average',5);
bgPix = filter2(fh,bgPix);
bgPix = bgPix > 0.5;

% Draw the figure
fig.fig = figure('Name','TQ-DESI',...
    'Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);

% Draw the ancilliary stuff
[fig] = drawAxes(fig,all,numC);

% Set the callback functions here for the individual things...
for n = 1:numC
    set(fig.sel(n,1),'Callback',{@singlePlot,fig,all,n,bgPix});
    singlePlot([],[],fig,all,n,bgPix);    
end

% Here we set the callbacks for the combi plot
for n = 1:3
    set(fig.op(n,1),'Callback',{@combiPlot,fig,all,bgPix});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig] = drawAxes(fig,all,qty)

% What is the size of the image?
sz = all(1).sz;
tmpImg = rand(sz(1),sz(2));

% Draw sub-axes down the left half of the figure
posn = fliplr(linspace(0.05,0.95,qty+2));

fig.ax = zeros(qty,2);
fig.sel = zeros(qty,1);
for n = 1:qty+1
    
    if n == qty + 1
        fig.main(1,1) = axes('Parent',fig.fig,...
            'Position',[0.05 posn(n+1)+0.01 0.6 posn(n)-posn(n+1)-0.01]);
        
        fig.main(1,2) = imagesc(tmpImg);
        
        set(fig.main(1),'XTick',[],...
            'YTick',[],...
            'XLim',[1 sz(2)],...    
            'YLim',[1 sz(1)],...
            'YDir','normal',...
            'LineWidth',5);
        box on;

        ylabel(fig.main(1),'Combined',...
            'FontSize',18,...
            'FontWeight','bold');

        
    else
        fig.ax(n,1) = axes('Parent',fig.fig,...
            'Position',[0.05 posn(n+1)+0.01 0.6 posn(n)-posn(n+1)-0.01]);
        
        fig.ax(n,2) = imagesc(tmpImg);
        
        set(fig.ax(n,1),'XTick',[],...
            'YTick',[],...
            'XLim',[1 sz(2)],...    
            'YLim',[1 sz(1)],...
            'YDir','normal',...
            'LineWidth',5);
        box on;

        ylabel(fig.ax(n,1),all(n).tqd.name,...
            'FontSize',18,...
            'FontWeight','bold');
        
        % Also draw me a listbox in which to display the various ions from
        % this image
        fig.sel(n,1) = uicontrol('Parent',fig.fig,...
            'Units','normalized',...
            'Position',[0.66 posn(n+1)+0.01 0.05 posn(n)-posn(n+1)-0.01],...
            'Style','listbox',...
            'String',num2str(all(n).tqd.head),...
            'FontSize',14,...
            'Min',1,...
            'Max',3,...
            'Value',1:numel(all(n).tqd.head));
        
    end
    
end

colormap(redbluecmap);


fig.bp1 = axes('Parent',fig.fig,...
    'Position',[0.75 0.6 0.245 0.35]);

fig.bp2 = axes('Parent',fig.fig,...
    'Position',[0.75 0.1 0.245 0.35],...
    'Visible','on');


% These for the dropdown boxes
vals = [0:1:qty]';
labs = cell(qty+1,1);
labs{1,1} = 'None';
for n = 1:qty
    labs{n+1,1} = all(n).tqd.name;
end

fig.op(1,1) = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.05 0.005 0.1 0.025],...
    'Style','popupmenu',...    
    'String',labs,...int2str(vals),...
    'Value',2);

if qty < 2
    tmp = 1;
else
    tmp = 3;
end
fig.op(2,1) = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.15 0.005 0.1 0.025],...
    'Style','popupmenu',...    
    'String',labs,...int2str(vals),...
    'Value',tmp);

if qty < 3
    tmp = 1;
else
    tmp = 4;
end
fig.op(3,1) = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.25 0.005 0.1 0.025],...
    'Style','popupmenu',...    
    'String',labs,...int2str(vals),...
    'Value',tmp);

fig.log = uicontrol('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0.9 0.01 0.1 0.05],...
    'Style','checkbox',...
    'String','Log?',...
    'Value',0);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function singlePlot(~,~,fig,all,i,bgPix)
% Draw a single plot in the corresponding place

tt = tic;

% Determine the ions that were selected...
ionInd = get(fig.sel(i,1),'Value');
%ionVal = str2num(get(fig.sel(i,1),'String'));
%ionVal = ionVal(ionInd);

% What is the correct image?
img = all(i).img(:,:,ionInd);

% Save the ion image as UserData of the listbox
set(fig.sel(i,1),'UserData',sum(img,3));

% Format for the single ion image display
img = sum(img,3);
img = img / max(img(:));

% Set all background values to 0
img(~bgPix) = 0;

% Take the image's sum and then replicate so that it is RGB whitescale
%img = repmat(img,[1 1 3]);

toc(tt);

% Plot it in the axes
set(fig.ax(i,2),'CData',img);

box on;

toc(tt)

% Now how about a mega box plot of all intensities included?
combiBoxPlot(fig,all);

% Update the actual main combi plot
try
    combiPlot([],[],fig,all,bgPix);
catch err
    err
end

toc(tt);

%drawnow;

%toc(tt);
disp('done');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function combiBoxPlot(fig,all)
% Draw a box plot of all selected ions - no scaling, just plain
clc;

numC = size(fig.ax,1);

globe = cell(numC,2);

% For each image, get the selected intensities...
for n = 1:numC
    
    % Determine the ions that were selected...
    ionInd = get(fig.sel(n,1),'Value');
    ionVal = str2num(get(fig.sel(n,1),'String'));
    ionVal = ionVal(ionInd);
    ionCll = cell(1,numel(ionInd));
    for r = 1:numel(ionCll)
        ionCll{1,r} = sprintf('%0.3f',ionVal(r));
    end

    % These are the selected ions
    img = all(n).img(:,:,ionInd);
    
    % Replicate the m/z vectors AND file name!
    sz = size(img);
    grp1 = repmat(ionCll,[sz(1)*sz(2) 1]);
    grp1 = grp1(:);    
    grp2 = repmat({all(n).tqd.name},[numel(grp1) 1]);
    
    % Reshape the bugger to make it a column
    %bpd = img(:);
    
    % Save all into global variable
    globe{n,1} = img(:);
    globe{n,2} = [grp2 grp1];

end

% Now let's try to make a boxplot
bpd = vertcat(globe{:,1});
lab = vertcat(globe{:,2});

bpd(bpd <= 0) = NaN;

axes(fig.bp1);
cla reset;
hold on;

boxplot(bpd,{lab(:,1) lab(:,2)},...
    'Colors','rbgm',...
    'Jitter',0.75,...
    'Symbol','o',...
    'LabelVerbosity','majorminor',...
    'FactorSeparator',[1],...
    'LabelOrientation','inline');%,...


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function individualPlot(ax,img,txt)
% Plot an individual image

axes(ax);
hold on;

imagesc(img);


colormap(gray);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function combiPlot(src,event,fig,all,bgPix)
% Draws the main plot...


% First get the callback values from the op boxes
tmp = cell2mat(get(fig.op,'Value'))-1;
labs = get(fig.op,'String');
labn = cell(numel(tmp),1);

% Now we need to create the correctly scaled image...
szI = size(all(1).img);
newImg = zeros(szI(1),szI(2),3);

dims = false(3,1);

for n = 1:3
    
    if tmp(n) ~= 0
        
        xx = get(fig.sel(tmp(n),1),'UserData');
        mx = max(xx(:));
        
        newImg(:,:,n) = xx / mx;
        
        dims(n,1) = true;
        
    end
    
    labn{n,1} = labs{n}{tmp(n)+1};
    
end

% Here we need to determine the background...

% This is the special method for Luisa - will override the existing
% controls...
% newImg2 = newImg;
% newImg2(:,:,2) = (newImg(:,:,1) + newImg(:,:,2)) / 2;
% newImg2(:,:,1) = 0;
% 
% newImg2(:,:,1) = newImg(:,:,3);
% newImg2(:,:,3) = 0;
% 
% % Divide green by red, put in green
% newImg2(:,:,2) = newImg2(:,:,2) ./ newImg2(:,:,1);
% tmp2 = newImg2(:,:,2);
% 
% newImg2(isinf(newImg2)) = 0;
% newImg2(isnan(newImg2)) = 0;
% 
% tmp2(tmp2 > 2) = 2;
% newImg2(:,:,2) = tmp2 ./ max(tmp2(:));
% 
% % Set all bg values to 0
% for n = 1:3
%     tmp = newImg2(:,:,n);
%     tmp(~bgPix) = 0;
%     newImg2(:,:,n) = tmp;
% end

% Create the labels
% labs = cell(szI(3)+1,1);
% labs{1,1} = 'N/A';
% for n = 2:szI(3)+1
%     labs{n,1} = sprintf('%0.3f',tqd.head(n-1));
% end

%newImg2(isnan(newImg2)) = 0;

% Which image to use?
%newImg = newImg2;
newImg = newImg;

% axes(fig.main);
% cla reset;
% hold on;
% %set(fig.main,'YDir','reverse');

% imagesc(flipud(newImg));

% set(gca,'XTick',[],'YTick',[],'YDir','reverse');

set(fig.main(2),'CData',newImg);

set(fig.main(1),'XTick',[],...
    'YTick',[],...
    'XLim',[1 size(newImg,2)],...    
    'YLim',[1 size(newImg,1)],...
    'YDir','normal',...
    'LineWidth',5);


%sz = size(newImg);
%xlim([0.5 sz(2)+0.5]);
%ylim([0.5 sz(1)+0.5]);

return

ylabel(fig.main(1),'Combined','FontSize',18,'FontWeight','bold');

% Histogram / box plot for the variance of the pixels...
axes(fig.bp2);
hold on;

%if sum(double(dims)) > 1
    set(fig.bp2,'Visible','on');
    [hh,subax,~,pp] = plotmatrix(reshape(newImg,[size(newImg,1)*size(newImg,2) size(newImg,3)]));
    
    set(hh,'Color','k','Marker','+');
    
    set(pp(1),'FaceColor','red');
    set(pp(2),'FaceColor','green');
    set(pp(3),'FaceColor','blue');
    
    ylabel(subax(1,1),labn(1),'FontWeight','bold','FontSize',12);
    ylabel(subax(2,1),labn(2),'FontWeight','bold','FontSize',12);
    ylabel(subax(3,1),labn(3),'FontWeight','bold','FontSize',12);
    
%else
    %set(fig.bp2,'Visible','off');
   
%end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

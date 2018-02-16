function makeConfMat(prc,qty,cols,order,txt)
% Make a uniform looking confusion matrix

% Calculate percentages if not provided
if isempty(prc)
    prc = 100 * bsxfun(@rdivide,qty,sum(qty,2));
end

% Size of the matrices
sz = size(prc,1);

if isempty(order)
    order = 1:sz;
end

% Text labels...
if isempty(txt)
    txtFlag = false;
    scSz = 200;
else
    txtFlag = true;
    txt = txt(order,:);
    scSz = 200;
end

% Let's allow us to pass cols as a 'groupInfo' cell, which contains colours
% (:,1) and symbols (:,2)... Might make life easier...
if iscell(cols)
    symb = vertcat(cols(:,2));
    cols = vertcat(cols{:,1});
else
    symb = [];
end

% We need to try to sort the order of the confusion matrices
prc2 = zeros(size(prc));
qty2 = zeros(size(qty));
for n = 1:sz
    for r = 1:sz
        prc2(n,r) = prc(order(n),order(r));
        qty2(n,r) = qty(order(n),order(r));
    end
end

if size(prc,2) == sz + 1
    prc2(:,end) = prc(order,end);
    qty2(:,end) = qty(order,end);
end

prc = prc2;
qty = qty2;

cmap = flipud(gray(100));

if sz >= 5
    fs = 18;
else
    fs = 24;
end

figure; hold on;

imagesc(prc);
colormap(cmap);
caxis([0 100]);

if ~txtFlag
    set(gca,'YDir','reverse','XTick',[],'YTick',[]);
else
    set(gca,'YDir','reverse',...
        'XTick',[],...
        'YTick',1:sz,...
        'YTickLabels',txt,...
        'FontSize',18);
    
end


if ~isempty(cols)
    
    % Symbols can be provided in this new version, but only need to bother
    % if we have provided them...
    x1 = 1:sz;
    y1 = zeros(sz,1) + 0.5;
    x2 = zeros(sz,1) + 0.5;
    y2 = 1:sz;
    
    if isempty(symb)
        
        % Just a single symbol, so easy
        scatter(x1,y1,scSz,cols(order,:),'o','filled',...
            'MarkerEdgeColor',[0.5 0.5 0.5]);
        scatter(x2,y2,scSz,cols(order,:),'o','filled',...
            'MarkerEdgeColor',[0.5 0.5 0.5]);

    else
        % Need to plot each individually...
        [ss,~,si] = unique(symb);
        for n = 1:numel(ss)
            fx = si == n;
            sy = si(fx);
        
            scatter(x1(fx),y1(fx),scSz,cols(order(fx),:),ss{n},...
                'filled','MarkerEdgeColor',[0.5 0.5 0.5]);

            scatter(x2(fx),y2(fx),scSz,cols(order(fx),:),ss{n},...
                'filled','MarkerEdgeColor',[0.5 0.5 0.5]);
        end
        
    end
end

for n = 1:sz
    
    for r = 1:sz
        
        % Skip empty pixels
        if qty(r,n) == 0
            continue;
        end
        
        
        if prc(r,n) > 60
            tc = 'white';
        else
            tc = 'black';
        end
        
        tl = [sprintf('%0.1f',prc(r,n)) '%' char(10) sprintf('%d', qty(r,n))];
        %tl = [sprintf('%0.1f',prc(r,n)) '%'];% char(10) sprintf('%d', qty(r,n))];
        
        text(n,r,tl,...
            'Color',tc,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...            %'FontUnits','pixels',...
            'FontSize',fs,...
            'FontWeight','bold');
        
    end
    
end

% Add unclassified pixels here
if size(prc,2) == sz + 1
    for r = 1:sz
        % Skip empty pixels
        if qty(r,end) == 0
            continue;
        end
                
        if prc(r,end) > 60
            tc = 'white';
        else
            tc = 'black';
        end
        
        tl = [sprintf('%0.1f',prc(r,end)) '%' char(10) sprintf('%d', qty(r,end))];
        
        text(sz+1,r,tl,...
            'Color',tc,...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...            %'FontUnits','pixels',...
            'FontSize',24,...
            'FontWeight','bold');
        
    end
end        
    
% Axis formatting stuffhere
ylabel('Actual Class   ','FontSize',20,'FontWeight','bold');
title('Predicted Class   ','FontSize',20,'FontWeight','bold');
box on;
axis tight square
set(gca,'LineWidth',5,...
    'TickLength',[0 0]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

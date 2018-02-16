function statsUnivariateStem(~,~,fig,window)
% Stats UV stem

% Get the guidata
sts = guidata(fig.fig);
if ~isfield(sts,'res') || ~isfield(sts.res,'uv')
    return
end

% Where are we drawing the things
parent = fig.ax.load;
parent2 = fig.ax.conf;

% What are the colours?
unqG = unique(sts.res.uv.grp);
numG = numel(unqG);
cols = parula(numG);

% Determine useful thresholds - eventually these will be defined by the
% user...
pqThresh = -log10(str2double(get(window.pqThresh,'String')));
fcThresh = str2double(get(window.fcThresh,'String'));

% Plot p or q values?
pqChoice = get(window.plotPQ,'Value');
if pqChoice == 1
    pqString = 'p';
else
    pqString = 'q';
end

% Log the pq values and use these throughout.  Need to remove the Inf
% values, or set them to be above the others by a slight way
pq = sts.res.uv.pq(:,pqChoice);
fc = sts.res.uv.fc;
lowPQ = min(pq(pq > 0));
idxLo = pq == 0;
pq(isnan(pq)) = 1;
pq = -log10(pq);
pq(idxLo) = -log10(lowPQ) * 1.2;

% Do below but for the plot above
f0 = get(parent2,'Children');
title(parent2,'');
delete(f0);

% Delete existing objects
f0 = get(parent,'Children');
delete(f0);

% Hold the axes
axes(parent);
hold on;

% Highlight only significant ones
fx = abs(fc) >= fcThresh & pq >= pqThresh;
%fx = pq > pqThresh;
if sum(fx) == 0
    fx = false(size(pq));
end


% Decide which values are to be plotted...
switch sts.datatype
    case 'ms'

        % Determine the axes limits and then draw the gridlines
        xl = [min(sts.proc.var.mz) max(sts.proc.var.mz)];
        yl = [-0.1 max(pq) * 1.05];

        % Horizontal line for PQ significance
        line([xl(1) xl(2)],[pqThresh pqThresh],...
            'LineStyle','--',...
            'Color','k');

        % Add in another for a higher value...
        line([xl(1) xl(2)],[floor(max(pq)) floor(max(pq))],...
            'LineStyle',':',...
            'Color','k');

        % Bad variables
        scatter(sts.proc.var.mz(~fx),pq(~fx),...
            60,[0.2627 0.5765 0.7647],'o','filled',...
            'MarkerEdgeColor',[0.7 0.7 0.7]);

        % Good variables
        scatter(sts.proc.var.mz(fx),pq(fx),...
            120,[0.8392 0.3765 0.3020],'o','filled',...
            'MarkerEdgeColor',[0.7 0.7 0.7],...
            'HitTest','on',...
            'ButtonDownFcn',{@statsUnivariateClick,...
            [parent parent2],...
            fc(fx),pq(fx),...
            sts.proc.sp(:,fx),...
            sts.res.uv.grp,...
            sts.proc.var.mz(fx),...
            cols,...
            pqString});
        
        % Axes info
        set(gca,...
            'YTickLabel',{[pqString ' = ' num2str(10^-pqThresh)],...
                [pqString ' = ' num2str(10^-floor(max(pq)))]},...
            'YTick',[pqThresh max([floor(max(pq)) pqThresh+1])],...
            'LineWidth',5,...
            'TickLength',[0 0],...
            'YDir','normal',...
            'XTick',[]);


    case 'lcms'
        
        % Determine the axes limits and then draw the gridlines
        xl = [min(sts.proc.var.mz) max(sts.proc.var.mz)];
        yl = [min(sts.proc.var.rt) max(sts.proc.var.rt)];

        % Bad variables
        scatter(sts.proc.var.mz(~fx),sts.proc.var.rt(~fx),...
            40,[0.2627 0.5765 0.7647],'o','filled',...
            'MarkerEdgeColor',[0.7 0.7 0.7]);
        
        % Good variables
        scatter(sts.proc.var.mz(fx),sts.proc.var.rt(fx),...%pq(fx),...
            120,[0.8392 0.3765 0.3020],'o','filled',...
            'MarkerEdgeColor',[0.7 0.7 0.7],...
            'HitTest','on',...
            'ButtonDownFcn',{@statsUnivariateClick,...
            [parent parent2],...
            fc(fx),pq(fx),...
            sts.proc.sp(:,fx),...
            sts.res.uv.grp,...
            [sts.proc.var.mz(fx) sts.proc.var.rt(fx) ...
            sts.proc.var.mz(fx) sts.proc.var.rt(fx)],...
            cols,...
            pqString});

        % Limited grid info
        set(gca,...
            'YTickLabel','',...
            'YTick',[],...
            'LineWidth',5,...
            'TickLength',[0 0],...
            'YDir','normal',...
            'XTick',[]);

end
    
box on;
grid off;

% Set the axes limits
xlim(xl);
ylim(yl);

% Link the bottom two x-axes
switch sts.datatype
    case 'ms'
        linkaxes([sts.fig.ax.spec parent],'x');
    case 'lcms'
        linkaxes([sts.fig.ax.spec parent],'xy');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

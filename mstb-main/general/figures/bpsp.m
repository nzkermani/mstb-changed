function [ output_args ] = bpsp(mz,var,grp,inds)
%bpsp - boxplotsubplot, provide a bunch of spectral variables, the grouping
%variable (e.g. histID) and a matrix of inds of how to lay out the subplots

[nr,nc] = size(inds)

figure; hold on;

yl = zeros(numel(inds),2);
n = 0;
for r = 1:nr    
    for c = 1:nc
        
        % Subplot counter
        n = n + 1;
        ax(n) = subplot(nr,nc,n);
        
        % This is the title
        tit = ['m/z = ' sprintf('%0.4f',mz(inds(r,c)))];
        
        % Actually draw it
        bp(ax(n),var(:,inds(r,c)),grp,tit);
        
        % Here determine the ylims for each on...
        yl(n,:) = [min(var(:,inds(r,c))) max(var(:,inds(r,c)))];
        
        %title([int2str(inds(r,c)) '-' int2str(n)]);
    end
end

set(ax,'YLim',[min(yl(:,1)) max(yl(:,2))]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bp(h,ydata,grp,tit)
% Draw the boxplot plot, calculate its mean and all that stuff...

% Determine x positions according to grouping variable
[unqN,~,unqI] = unique(grp);

% Add a little wiggle to the data...
wigval = 0.5;
wiggle = rand(numel(unqI),1) * wigval;
wiggle = wiggle - wigval/2;
xval = unqI + wiggle;

% Scatter the points
scatter(xval,ydata,50,'blue','o');

% Add the mean and median values for each group...
for n = 1:numel(unqN)
    
    % Indices
    fx = find(unqI == n);
    
    % Mean
    mn = mean(ydata(fx,:));
    
    % Median
    md = median(ydata(fx,:));
    
    % Draw on plot...
    line([n-wigval/2 n+wigval/2],[mn mn],...
        'LineWidth',4,'Color','red');
    
    line([n-wigval/2 n+wigval/2],[md md],...
        'LineWidth',4,'Color','magenta');
end

title(tit,'FontSize',18);
xlim([min(unqI)-0.5 max(unqI)+0.5]);
set(gca,'XTick',1:numel(unqN),'XTickLabel',unqN)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% % Use this data...
% ydata = op.XPeaksLog(:,px);
% 
% 
% 
% legend(hnd,{'Mean','Median'},'Location','NorthEastOutside',...
%     'Orientation','vertical');
% 
% ylabel('Logged intensity','FontSize',14);
% 


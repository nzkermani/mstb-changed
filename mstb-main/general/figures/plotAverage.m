function plotAverage(x,y,grps,cols)

% Unique groups
[unq,~,ind] = unique(grps);
numG = numel(unq);
h = zeros(numG,1);

% Correct for numeric grps
if isnumeric(unq)
    nunq = cell(size(unq));
    for n = 1:numG
        nunq{n,1} = num2str(unq(n));
    end
    unq = nunq;
end

% Create new figure
figure('Position',[200 200 1000 500]); hold on;

% Loop through
for n = 1:numG
    
    % Which obs?
    fx = ind == n;
    
    % Determine percentiles
    prc = prctile(y(fx,:),[25 50 75]);    

    % These are the boundaries of the 1st 3rd
    vals = [prc(1,:)' prc(3,:)'-prc(1,:)'];

    % Plot the boundary
    hn = area(x,vals);
    
    % Remove the excess one
    set(hn(1),...
        'FaceColor','none',...
        'EdgeColor','none');
    
    % Change the colour of the new one
    set(hn(2),...
        'FaceColor',cols(n,:),...
        'EdgeColor','none');
    
    % Set the transparency
    %set(get(hn(2),'Children'),'FaceAlpha',0.1);
    hn(2).FaceAlpha = 0.25;
    
    % Plot the line
    plot(x,prc(2,:),...
        'LineWidth',2,...
        'Color',cols(n,:));
    
    h(n,1) = hn(2);

end

xlim([min(x) max(x)]);

% Legend
legend(h,unq);

set(gca,'FontSize',16);
xlabel('Wavenumber / cm-1','FontSize',18,'FontWeight','bold');
ylabel('Intensity','FontSize',18,'FontWeight','bold');
box on;


end
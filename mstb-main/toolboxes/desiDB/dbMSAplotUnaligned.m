function dbMSAplotUnaligned(data)
%dbMSAplotUnaligned - plot the average of each tissue type from each file.

numF = size(data,2);

% Somewhere to store everything, but which could expand if necessary
avgs = cell(100,3);
i = 0;
fileNo = zeros(100,1);

for n = 1:numF
    
    [hist,~,idx] = unique(data(n).histID);
    numG = numel(hist);
    
    for r = 1:numG
        
        % Indices
        fx = idx == r;
       
        % Calculate the mean spectrum
        i = i + 1;
        avgs{i,3} = nanmean(data(n).sp(fx,:),1);
        avgs{i,2} = data(n).mz;
        
        % Label for legend
        try
            avgs{i,1} = [data(n).file(1:end-4) '-' hist{r}];
        catch
            avgs{i,1} = [data(n).name '-' hist{r}];
        end
                
        
        % File number
        fileNo(i,1) = n;
        
    end
    
end

% Ditch extra entries in the cell
avgs = avgs(1:i,:);

% Calculate colours
cols = flipud(parula(i));




% Draw the figure
fig = figure('Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5]);

hold on;
ax = zeros(i,1);
for n = 1:i
    
    [a,b] = insertZeros(avgs{n,2},avgs{n,3},0.005);
    b = 1e4 * b / sum(b);
    
    b = b + ( (fileNo(n,1)-1) * 100 );
    
    ax(n,1) = plot(a,b,...
        'Color',cols(n,:),'LineWidth',2);
    
end

legend(ax,avgs(:,1))

    

end


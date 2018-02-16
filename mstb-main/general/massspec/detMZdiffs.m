function [ fit ] = detMZdiffs(data,fig)
%detPPMdiffs - the function that will scan through the orbData strcuture
%thing and calculate Th/Da differences between neighbouring data points.

if nargin == 1
    fig = false;
end

% Define a Th/Da threshold at which point we decide points are not contiguous
thresh = 0.1;

% Structure in which to store the results
res = struct('data',[],'inds',[]);
c = 0;

% How many observations in the dataset?
numO = size(data,2);
if size(data,1) > numO
    numO = size(data,1);
    flag = true;
else
    flag = false;
end

% What kind of fitting to do. Quadratic required at least
fitVal = 1;

% If we are plotting things...
if fig
    tmpFig = figure('Units','normalized','Position',[0.25 0.25 0.5 0.5]);
    ax(1) = gca;%subplot(1,2,1); hold on;
    %ax(2) = subplot(1,2,2); hold on;
    hold on;
    %cols = jet(numO);
end


% Loop through each 'observation'
for n = 1:numO
    
    % How many scans in the observation?
    if flag
        numS = 1;
    else
        numS = size(data{n},2);
    end
    
    % Save the coefficients of individual scans
    if fig
        plotInd = sort(randperm(numS,min([numS 3])));        
    end
    
    % Loop through each scan
    for s = 1:numS
                
        % Get the mz vector
        if flag
            mzv = data{n}(:,1);
            inv = data{n}(:,2);
        else
            mzv = data{n}{s}(:,1);
            inv = data{n}{s}(:,2);
        end
        
        % Determine the sampling frequency over points        
        mzdiff = [diff(mzv); NaN];
        
        % Convert to ppm values
        %mzdiff = mzdiff ./ mzv .* 1e6;
        
        % Find diffs that are less than the threshold
        fx = mzdiff <= thresh & mzdiff > 0 & inv > 0;
        
        c = c + 1;
        res(c).data = [mzv(fx) mzdiff(fx)];
        res(c).inds = [n s];
        
        % Do individual fitting here to see if there are differences...
%         if fig && any(plotInd == s)           
%             
%             % Individual fit for this scan
%             loopFit = polyfit(mzv(fx),mzdiff(fx),fitVal);
%             
%             % New x and y values
%             tmpX = min(mzv(fx)):1:max(mzv(fx));
%             tmpY = polyval(loopFit,tmpX);
%             
%             % Add to the figure
%             plot(ax(1),tmpX,tmpY,'--','Color',cols(n,:));
%         end 
    end
        
end
        
% Can we concatenate the results into a single matrix...
all = vertcat(res.data);

% Sort
all = sortrows(all,1);

% Because the fitting is dominated by more frequently ocurring peaks, we
% should ensure to sample peaks periodically so that we don't favour m/z
% regions with more peaks, but rather favour all regions equally. So this
% requires some ditching of m/z values that occur too frequently...
res = 0.5;
rndX = round(all(:,1) ./ res);
rndZ = (rndX(1) * res):res:(rndX(end) * res);
rndX = rndX - rndX(1) + 1;
rndY = accumarray(rndX,all(:,2)',[],@median);

idx = rndY > 0;

% Now we have the regular function that can be fit, without being overfit
% by virtue of density of datapoints
maxPoly = 1;
cols = jet(maxPoly);
if maxPoly == 3
    cols = [0 0 0; 1 0 0; 0 0 1];
end

fit = struct('p',[],'s',[],'mu',[]);
hand = zeros(maxPoly,1);
labs = cell(1,maxPoly);

for n = 1:maxPoly
    
    %[fit(n).p,fit(n).s,fit(n).mu] = polyfit(rndZ',rndY,n);
    [fit(n).p] = polyfit(rndZ(idx)',rndY(idx),n);
    
    
    [newy] = polyval(fit(n).p,rndZ);%,fit(n).s,fit(n).mu);
    
    if fig
        hand(n,1) = plot(ax(1),rndZ,newy,'Color',cols(n,:),'LineWidth',10);
        labs{n} = ['  Order = ' int2str(n)];
    end
end

if fig
    idx2 = randperm(numel(idx),numel(idx)-min([numel(idx) 200]));
    idx(idx2) = false;
    hand(1,1) = scatter(rndZ(idx),rndY(idx),20,'w','o','filled',...
        'MarkerEdgeColor','k');
    labs{1} = '  Raw subset';
    xlabel('m/z','FontSize',16,'FontWeight','bold','FontAngle','italic');
    ylabel('Difference / Th','FontSize',16,'FontWeight','bold');
    legend(hand,labs,'Location','NorthWest');
    box on;
    set(gca,'FontSize',14);
    
    figure('Units','normalized','Position',[0.25 0.25 0.5 0.5]); hold on;
    idx = randperm(size(all,1),min([5000 size(all,1)]));
    scatter(all(idx,1),all(idx,2),30,'o','filled');
    xlabel('m/z','FontSize',16,'FontWeight','bold','FontAngle','italic');
    ylabel('Difference / ppm','FontSize',16,'FontWeight','bold');
    box on;
    set(gca,'FontSize',14);

end



% This is what to return
fit = fit(fitVal).p;


end


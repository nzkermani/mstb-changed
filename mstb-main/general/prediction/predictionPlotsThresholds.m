function [tab] = predictionPlotsThresholds(data)
%imagePlots - all sorts of plots for Jocelyn

threshs = 0.1:0.1:1;

[tab] = gridPlot(data,threshs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tab] = gridPlot(data,threshs)

images = [1 2 4];
bg = 3;

fig.fig = figure('Units','normalized',...
    'Position',[0.25 0.0 0.2285 0.8868]);

numT = numel(threshs) + 1;
numD = size(data,2);
ax = zeros(numT,numD);

% Draw the axes
for n = 1:numT
    
    % Axes height
    ht = (numT - n) / numT;
    
    for r = 1:numD
        
        % Axes left position
        lt = (r-1) / numD;
        
        ax(n,r) = axes('Parent',fig.fig,...
            'Units','normalized',...
            'Position',[lt ht 1/numD 1/numT]);
        
    end
    
end

% Somewhere to store the pixel classification rates
tab = struct('file',[],'correct',[],'quantity',[],'labels',[]);

% Now for the subplots...
for n = 1:numD
    
    % Get the image
    try
        tmpUntouch = data(n).prob(:,:,images);
    catch
        tmpUntouch = data(n).Pr(:,:,images);
    end
    
    % Reshape TOBG if necessary...
    if size(data(n).tobg,2) == 1
        tobg = reshape(data(n).tobg,data(n).sz(1:2));
    else
        tobg = data(n).tobg;
    end    
    bgprob = data(n).Pr(:,:,bg);    
    tobg = (tobg + (1-bgprob)) == 2;
       
    
    % What about the annotations for each of the images at each of the
    % thresholds.  Need to extract the annotations from each file
    pr2 = reshape(data(n).Pr,[data(n).sz(1)*data(n).sz(2) size(data(n).Pr,3)]);
    annoIdx = data(n).anno.mask > 0;
    pr2 = pr2(annoIdx,:);
    
    % Determine the number of annotated pixels in this file (for each
    % class)
    [classID,~,classInd] = unique(data(n).anno.histID(annoIdx));
    classQty = zeros(size(pr2,2),1);
    for c = 1:size(pr2,2)
        classQty(c,1) = sum(classInd == c);
    end
        
    % Huge confusion matrix of results
    cmat = zeros(numel(classID),numel(classID)+1,numT-1);
    
    % Do differently for each threshold
    for t = 1:numT
        
        % Determine the confusion matrix
        if t < numT            
            [cmat(:,:,t)] = confusionMatrix(classInd,pr2,threshs(t));
        end
        
        % Get the 'pure' prob image
        tmp = tmpUntouch;
        
        % Apply the TOBG mask to the image
        for r = 1:3
            
            % This is the channels pr
            ff = tmp(:,:,r);
            
            % Threshold it!
            if t < numT
                ff = double(ff >= threshs(t));
            else
                ff = double(ff);
            end
            
            % Set the background to another colour
            ff(~tobg) = 0.9;
            
            % Save back to the image
            tmp(:,:,r) = ff;
        end
        
        % Set the axes
        axes(ax(t,n));
        
        imagesc(tmp);
        
        set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none');
        
        % Add text label
        if n == 1 && t < numT
            text(5,15,sprintf('%0.1f',threshs(t)),...
                'FontSize',16,...
                'FontWeight','bold',...
                'Color','black');
            
        elseif n == 1 && t == numT
            text(5,15,'0-1',...
                'FontSize',16,...
                'FontWeight','bold',...
                'Color','black');
            
        end
    end
    
    % Need to do something with the pixel classification rates
    tab(n).quantity = classQty';
    tab(n).labels = classID;
    tab(n).file = n;
    tab(n).cmat = cmat;
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

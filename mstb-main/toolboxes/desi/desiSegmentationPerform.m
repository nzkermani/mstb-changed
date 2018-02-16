function desiSegmentationPerform(~,~,fig,seg)
%desiSegmentationPerform - calculate the segmentatation requested and work
%out how to show it in the figures...

% Get the data
dpn = guidata(fig.fig);

% Now determine which of the methods is to be performed
methods = get(seg.list,'String');
val = get(seg.list,'Value');
choice = methods{val};

% Ignore the background pixels?
bg = get(seg.remBG,'Value');
if bg
    tobg = dpn.d1.tobg;
else
    tobg = [];
end

% We need to fill in the background image quite a lot in order to make it
% one homogeneous area rather than potentially having a lot of holes in it
if ~isempty(tobg)
    
    % Just define the options here - may need to make these interactive
    % when we move towards different sections
    mOpts.imopen = 2;
    mOpts.imclose = 8;
    mOpts.bigoper = [];
    
    % This does the operations defined.
    tobg = doMorph(tobg,mOpts);
end


% Do we apply a log transform to the data?
doLog = get(seg.doLog,'Value');

% Also determine the m/z values to be used
mzL = str2double(get(seg.mzL,'String'));
mzH = str2double(get(seg.mzH,'String'));

% Normalisation
methods = get(seg.norm,'String');
val = get(seg.norm,'Value');
norm = methods{val};

% Max clusters?
maxK = str2double(get(seg.maxClust,'String'));
if ~isnumeric(maxK)
    maxK = 4;
end

% Now we can just proceed to run the various functions...
switch choice    
        
    case {'PCA','kMeans','MMC'}
        
        % Extract the data matrix        
        [dx1,dx2,anno] = dpnPrepFullMS(dpn,[mzL mzH],tobg);
        
        
    otherwise
        
        % Return the images to the ion images...
        set(dpn.fig.ax.ms1(2),'CData',dpn.d1.img);
        set(dpn.fig.ax.ms1(1),...
            'XLim',[1 size(dpn.d1.img,2)],...
            'YLim',[1 size(dpn.d1.img,1)],...
            'XTick',[],...
            'YTick',[]);
                
        % Change the text...
        set(fig.txtPos,'String','Ion image');

        % Set the bg rem field to 0 to make consistent
        set(seg.remBG,'Value',0);
        
        % Remove functionality from the back/forward buttons
        set(seg.imgPrev,'Callback',[]);

        % Do this for the second axes if applicable
        if isfield(dpn,'ms2')
            set(dpn.fig.ax.ms2(2),'CData',dpn.d2.img);
            set(dpn.fig.ax.ms2(1),...
                'XLim',[1 size(dpn.d1.img,2)],...
                'YLim',[1 size(dpn.d1.img,1)],...
                'XTick',[],...
                'YTick',[]);        
            set(fig.txtNeg,'String','Ion image');
            set(seg.doLog,'Value',0);        
            set(seg.imgNext,'Callback',[]);
        end
        
        % Job done - return satisfied
        return
end
   
% Create a waitbar.
wb = waitbar(0,'Performing...');

% Perform the normalisation method as necessary...
switch lower(norm)
    
    case 'pqn-median'
        dx1.sp = jsmNormalise(dx1.sp,'pqn-median',0,[]);
        
    case 'pqn-mean'
        dx1.sp = jsmNormalise(dx1.sp,'pqn-mean',0,[]);
        
    case 'tic'
        ss1 = mean(nansum(dx1.sp,2));
        dx1.sp = jsmNormalise(dx1.sp,'tic',0,[]);        
        dx1.sp = dx1.sp * ss1;
        
    otherwise
        % No normalisation is to be performed
        
end

% If we are going to log the data, we need to determine the optimal offset
% to be applied to each dataset...
if doLog    
    % Offsets
    os1 = nanmedian(dx1.sp(dx1.sp > 0));
    
    % Transform
    dx1.sp = log(dx1.sp + os1);    
end

% Here is where we can be specific about the analysis to be performed. We
% have extracted the variables to be used above, so here is where we
% diversify
switch choice
    case 'PCA'
        
        % Generic waitbar progress
        waitbar(1/2,wb,'PCA');
                
        % Do PCA...
        [dx1.load,dx1.data,dx1.ee] = princomp(dx1.sp,'econ');
        dx1.ee = 100 * dx1.ee / sum(dx1.ee);
        
        dx1.data = dx1.data(:,1:10);
        
        dx1.l = cell(10,1);
        for n = 1:10
            dx1.l(n,1) = {['PC' int2str(n) ', ' sprintf('%0.2f',dx1.ee(n)) '%']};
        end
        
    case 'kMeans'
        
       
        % Generate empty matrices
        dx1.data = zeros(size(dx1.sp,1),maxK-1);
        dx1.l = cell(maxK-1,1);
        dx1.load = [];
                
        for n = 2:maxK
            
            % Waitbar
            waitbar(n/maxK,wb,'k-means');
            
            % Run the kmeans functions
            dx1.data(:,n-1) = kmeans(dx1.sp,n);
            
            % Write the little label
            dx1.l(n-1,1) = {['k = ' int2str(n)]};                       
        end
        
    case 'MMC'
        
        % Indices of the training set
        %[mask2,histID,isInt1,isInt2,pixID] = dpnAnnotationExtract(dpn);
        idx = anno.mask2 > 0;
        
        % What is the mean spectrum?
        mean1 = nanmean(dx1.sp(idx,:),1);
        
        % Subtract the mean.
        tmp1 = bsxfun(@minus,dx1.sp(idx,:),mean1);
       
        % Run the OAA-MMC function
        [B1,dx1.load] = oaaMMC(tmp1,anno.histID(idx),2);
        
        % Calculate probabilities for pixels
        [~,Pr1,~] = calcScoresSpectra(dx1.sp,mean1,dx1.load,B1);

        dx1.l = cell(size(Pr1,2)+1,1);
        for r = 1:size(Pr1,2)
            dx1.l(r,1) = {['LV ' int2str(r)]};
        end
        dx1.l(r+1,1) = {'Composite'};
            
        dx1.data = Pr1;
        
    otherwise
        % There is no otherwise
end

% So these results need to be saved somewhere in the guidata
[dpn.d1.mva.(choice).img] = return2image(real(dx1.data),dx1.idx,size(dpn.d1.sp));
dpn.d1.mva.(choice).label = dx1.l;
dpn.d1.mva.(choice).load = dx1.load;
dpn.d1.mva.(choice).doLog = doLog;
dpn.d1.mva.(choice).mz = [mzL mzH];
dpn.d1.mva.(choice).norm = norm;


% Save to the guidata
guidata(fig.fig,dpn);

% Here we want to update the boxes in the panel with the number of
% components (etc) that were calculated, and a 'NA' box at the top to
% ignore the channel
sz = size(dpn.d1.mva.(choice).img,3);
switch choice
    case 'kMeans'
        labs = cell(sz,1);
        for n = 1:sz
            labs{n,1} = int2str(n+1);
        end
        
        % Update the boxes in the GUI        
        set(seg.res,'String',labs);
        set(seg.res(1),'Value',1,'Visible','on');
        set(seg.res(2),'Value',1,'Visible','off');
        set(seg.res(3),'Value',1,'Visible','off');

    otherwise        
        % This is for everything else which has composite images
        labs = cell(sz+1,1);
        labs{1,1} = 'N/A';
        for n = 1:sz
            labs{n+1,1} = int2str(n);
        end
        
        % Update the boxes in the GUI
        set(seg.res,'String',labs);
        set(seg.res(1),'Value',min([2 numel(labs)]),'Visible','on');
        set(seg.res(2),'Value',min([3 numel(labs)]),'Visible','on');
        set(seg.res(3),'Value',min([4 numel(labs)]),'Visible','on');

end


% Now look to refresh the images on display...
%dpnIonImage([],[],dpn.fig.ax.mv,dpn.d1.mva.(choice).img(:,:,1:3));

% Change the text...
set(dpn.fig.txtPos,'String','PCs 1-3');%dpn.d1.mva.(choice).label{1});

% Set the callback function for the three RGB boxes
set(seg.res,'Callback',{@dpnChangeSegImage2,fig,seg,choice});
dpnChangeSegImage2([],[],fig,seg,choice);

% Delete the waitbar
delete(wb);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmpData] = return2image(data,index,sz)

% What is the size of the new images?
newSz = [sz(1) sz(2) size(data,2)];

% Full size empty matrix
tmpData = NaN([newSz(1)*newSz(2) newSz(3)]);

% Restore to original location
tmpData(index,:) = data;

% Reshape
tmpData = reshape(tmpData,newSz);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx1,dx2] = pcaImageCheck(dx1,dx2,choice)

numI = size(dx1.mva.(choice).img,3);

for n = 1:numI
    
    % Check addition and subtraction of the images
    a1 = dx1.mva.(choice).img(:,:,n) + dx2.mva.(choice).img(:,:,n);

    s1 = dx1.mva.(choice).img(:,:,n) - dx2.mva.(choice).img(:,:,n);
    
    hst = hist([a1(:) s1(:)],[-5:1:5]);
    
    % Decide which is the best based on the number of near-0 intensity
    % pixels
    if hst(6,1) > hst(6,2)
        dx2.mva.(choice).img(:,:,n) = dx2.mva.(choice).img(:,:,n) * -1;
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx1,dx2] = kMeansImageCheck(dx1,dx2,choice)
% Similar idea to the PCA function above, but needs to look at the clusters
% individually

% Short cuts to the images
im1 = dx1.mva.(choice).img;
im2 = dx2.mva.(choice).img;

% How many images to do this for...
numI = size(im1,3);
    
% Loop through each image
for n = 1:numI
    
    % How many levels (clusters) are there in this image?
    numC = max(max(im1(:,:,n)));
    
    crs  = zeros(numC,numC);
    % So for each cluster
    for r = 1:numC
        
        % Extract the main image
        x = double(im1(:,:,n) == r);
        y = im2(:,:,n);
        
        
        for p = 1:numC
            
            yy = double(y == p);
            crs(p,r) = corr(x(:),yy(:));
        end
    end
    
    % So now with each cluster we find the maximum value
    maxVal = max(crs,[],1);
    [~,idx] = sort(maxVal,'descend');
        
    % What values should be changed?
    crs2 = crs;
    map = zeros(numC,2);
    for r = 1:numC
        
        % Find the row of the largest value in column idx(r)
        [a,b] = max(crs(:,idx(r)));
        
        % Apply to the map
        map(r,:) = [b idx(r)];
        
        crs(b,:) = NaN;
                
    end
    
    % Now apply the map to the image...
    new2 = im2(:,:,n);
    for r = 1:numC
        fx = im2(:,:,n) == map(r,1);
        new2(fx) = map(r,2);
    end
    
    im2(:,:,n) = new2;

end

dx2.mva.(choice).img = im2;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

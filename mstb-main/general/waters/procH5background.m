function [final] = procH5background(mz,sp,sz)
%procH5background - determine the background for an individual section. We
%need to provide only the pixels of one image and its mz vector

% There are various ways to do this - it could be the otsu method, the
% kmeans approach.  Need to be able to decide which is best...
mzRange = 'non-raff';
doLog = true;
method = 'kmeans';

% Function handle for reshaping...
mat2img = @(x) reshape(x,[sz(1) sz(2) size(x,2)]);
img2mat = @(x) reshape(x,[sz(1)*sz(2) size(x,3)]);

% With waters data we have problems with the first few columns being either
% blank or at a phenomenal intensity
sp = mat2img(sp);
sp(:,1,:) = sp(:,2,:);
sp = img2mat(sp);

% Define the mz ranges for the various options...
switch mzRange
    case 'full'
        mask = mz > 0;        
    case 'non-raff'
        mask = ~(mz > 500 & mz < 600);        
    case 'lipid'
        mask = mz > 600;        
    otherwise
        % No otherwise
end

% Do we want to log the data?
if doLog
    os = nanmedian(sp(sp > 0));
    sp = log(sp + os);
end

% Now run the method...
switch method
    
    case 'otsu'
        
        tmp = nansum(sp(:,mask),2);        
        tobg = dpnTOBG(tmp,[],[]);
        
    case 'kmeans'
        
        tobg = kmeans(sp(:,mask),2,...
            'Distance','cosine',...
            'Replicates',3,...
            'OnlinePhase','on',...
            'MaxIter',200,...
            'Display','iter');
    case 'pca'
        
        [ll,ss,ee] = pca(sp(:,mask),...
            'NumComponents',3);
        
        tobg = dpnTOBG(ss(:,1),[],[]);        
        
    otherwise
        disp('No otherwise');
end

% Reshape to be an image
tobg = mat2img(tobg);

% If doing kmeans, need to determine the index of the tissue and set tissue
% to equal 1 and background 0
if strcmp(method,'kmeans')
    
    % The corners should be the best places to look...
    corners = false(size(tobg));
    corners(1:5,1:5) = true;
    corners(1:5,end-4:end) = true;
    corners(end-4:end,1:5) = true;
    corners(end-4:end,end-4:end) = true;
    
    % Determine the mode value
    modv = mode(tobg(corners));    
    tobg = tobg ~= modv;
    
elseif strcmp(method,'otsu')
    tobg = tobg == 1;    
end

% Smooth
smtobg = smoothImage(tobg,5);

% Determine the final threshold for the background
final = smtobg > 0.5;
final = img2mat(final);

return

% Plot...
figure('Units','normalized','Position',[0.1 0.1 0.4 0.6]); 
ax(1) = subplot(2,2,1); imagesc(mat2img(nansum(sp,2)));
ax(2) = subplot(2,2,2); imagesc(tobg);
ax(3) = subplot(2,2,3); imagesc(smtobg);
ax(4) = subplot(2,2,4); imagesc(smtobg > 0.5);
linkaxes(ax,'xy');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = smoothImage(x,sz)
% Simple image smoothing function

filt = fspecial('disk',sz);
y = filter2(filt,x);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

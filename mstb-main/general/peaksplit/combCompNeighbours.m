function [ mz,x ] = combCompNeighbours(mz,x,ppmTol,flag)
%combCompNeighbours - combine neighbouring images that are within a certain
%ppm tolerance, and have PERFECTLY complementary binary images.

verbose = true;

% In some instances we don't want the data matrices to change sizes, so we
% need to return the matrices in their original size. This can be
% accomplished by providing a fourth 'flag' input as TRUE
if nargin == 3
    flag = false;
end

% It would be good to add an extra empty image at top
% tmp = zeros(size(x,1),size(x,2));
% x = cat(3,tmp,x);
% try
%     mz = [0.1 mz];
% catch
%     mz = [0.1; mz];
% end

% Create a figure to show off the results
if verbose
    fig.fig = figure;
    fig.ax(1) = subplot(1,3,1); 
    fig.im(1) = imagesc(rand(size(x,1),size(x,2)));
    fig.ax(2) = subplot(1,3,2); 
    fig.im(2) = imagesc(rand(size(x,1),size(x,2)));
    fig.ax(3) = subplot(1,3,3); 
    fig.im(3) = imagesc(rand(size(x,1),size(x,2)));
    linkaxes(fig.ax,'xy');
end

% Number of variables
numV = size(x,3);

% List of the variables that are to be deleted
del = zeros(numV,2);

% Determine the mean of each variable...
tic
allMean = nansum(x,1);
allMean = nansum(allMean,2);
del(:,2) = squeeze(allMean);
%allMean = nanmean(reshape(x,[size(x,1)*size(x,2) size(x,3)]),1);
clear allMean
toc

for n = 1:numV-1
    
    % Determine if neighbours are complementary
    neigh = n + 1;
    
    % Now with the neighbour, determine the ppm difference
    ppmDiff = 1e6 * (mz(neigh) - mz(n)) / mz(n);
    %ppmChk  = abs(ppmDiff) <= ppmTol;
           
    if abs(ppmDiff) <= ppmTol
        
        
        % Add up the binary images
        imgs = double(x(:,:,n) > 0) + double(x(:,:,neigh) > 0);

        if verbose
            set(fig.im(1),'CData',x(:,:,n) > 0);
            set(fig.im(2),'CData',x(:,:,neigh) > 0);
            set(fig.im(3),'CData',imgs);
            drawnow;
        end
        
        % If the biggest value is 1, i.e. full complementarity/no overlap
        % then we can merge the variable into the neighbouring variable
        if max(imgs(:)) == 1
                        
            % Mark which variables are to be merged with 'del'
            del(n,1) = 1;
            %tmp = x(:,:,n);        
            
            % Combine the images
            x(:,:,neigh) = x(:,:,neigh) + x(:,:,n);
            x(:,:,n) = 0;
                        
            % Print to screen / waitbar would be better
            disp([int2str(n) ' - ' ...
                sprintf('%0.3f',mz(n)) ' - ' ...
                sprintf('%0.3f',mz(neigh))]);
                        
        end
        
        
        
    end
    
    
end

% Now we need to calculate the mean m/z value again
for n = 2:numV
    
    % Is this peak 0 and the one(s) before it deleted?
    if del(n,1) == 0 && del(n-1,1) == 1
                
        % Then find all previous values of 1, i.e. these are the previous
        % neighbours that were merged into this variable
        fx = find(del(1:n-1,1) == 0,1,'last') + 1;
        
        % So now from fx ... n are the variables that have been combined
        % together
        %mzvals = mz(fx:n);
        
        % This is the weighted m/z value representative of the intensities
        % of its constituent parts
        wtm = sum(mz(fx:n) .* del(fx:n,2)') / sum(del(fx:n,2));
        mz(n) = wtm;
        mz(fx:n-1) = NaN;        
        
    end
    
end

if verbose
    close(fig.fig);
end

% Set the first dummy variable to NaN so that we ensure it is deleted
%mz(1) = NaN;
   
% Trim out the variables that were deleted
fx = isnan(mz);

% Decide on how to present the output
if flag
    
    % Just for info
    disp(['>>> Variables removed = ' int2str(sum(fx))]);
    
    dummy = (1:sum(fx))/1000;
    mz(fx) = dummy;
    
else
    
    % Here we delete the variables
    mz = mz(~fx);
    x  = x(:,:,~fx);
    
end



end


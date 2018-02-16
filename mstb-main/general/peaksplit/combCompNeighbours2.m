function [ mz,x ] = combCompNeighbours2(mz,x,ppmTol,flag)
%combCompNeighbours2 - this is a new version which is supposed to be better
%than the first version.  It does this by merging the most likely ones
%before others; i.e. work on highly populated m/z bins first, and ignore
%ones which look most like noise.

percThreshold = 1;

% In some instances we don't want the data matrices to change sizes, so we
% need to return the matrices in their original size. This can be
% accomplished by providing a fourth 'flag' input as TRUE
if nargin == 3
    flag = false;
end

% For each variable, we need to determine the population of the image, i.e.
% how many non-zero intensities it has.  Also the mean intensity of
% non-zero pixels too.
numV = size(x,3);
info = zeros(numV,3);
for n = 1:numV
    
    fx = x(:,:,n);
    fx = fx(:);
    
    info(n,:) = [n sum(fx > 0) mean(fx(fx > 0))];
    
end
% figure; scatter(info(:,2),info(:,3),70,info(:,1));

% Now with that information we can rank the  variables which are better and
% which should be merged first... Let's sort on population
[srt,idx] = sortrows(info,-2);

% So could just work down the list of 'srt' to determine which variables
% are potentially merge-able with each of the most populous.  Start with
% one variable either side (within ppm tolerance) and extend only if
% complementary.  If n-1 and n+1 don't work together, then pick the one
% with the best fit.  Don't forget to set merged variables to zero to
% prevent them being used twice.
isMerged = false(numV,1);

for n = 1:numV
    
    % Always refer to its true index, rather than the sorted one
    i = srt(n,1);
    
    % Check that the variable hasn't been merged elsewhere
    if ~isMerged(i,1)
        
        % Mark this one as having been seen to...
        isMerged(i,1) = true;
        
        r = 0;
        ppmDiff = [0 0];
        flag = true;        
                       
        % Need to look at its nearest neighbours...
        while abs(min(ppmDiff)) < ppmTol && flag
            
            % This will tell us which to merge in the end
            doMerge = [false false];
            
            % Increase r
            r = r + 1;
            
            % Determine neighbours
            neigh = [-r r] + i;
            
            % Check that these variables actually exist
            vchk = neigh > 0 & neigh <= numV;
            neigh = neigh(vchk);
            
            % Have these been merged before?
            mergeCheck = isMerged(neigh,1)';
            
            % Determine their ppm difference
            ppmDiff = 1e6 * (mz(neigh)-mz(i)) / mz(i);
            
            % If the variables have been merged already, then spike their
            % ppm differences to be too large
            ppmDiff(mergeCheck) = Inf;
            
            % Check which fall within the tolerance
            ppmCheck = abs(ppmDiff) <= ppmTol;
            
            % So it may be the case that either/both/none of the adjacent
            % variables are within the tolerance
            if sum(ppmCheck) == 0
                % Then neiher are within, so skip
                flag = false;
                                
            elseif sum(ppmCheck) == 1
                % Just one could be included. Check that it is
                % complementary
                poss = neigh(ppmCheck == 1);
                
                img1 = double(x(:,:,i) > 0);
                img2 = double(x(:,:,poss) > 0);
                
                % Sum them together and see if the maximum value is 1
                imgX = img1 + img2;
                multiCount = sum(imgX(imgX > 1));
                multiPerc = 100 * multiCount / numel(imgX);                
                
                if multiPerc < percThreshold
                    % Then images are complementary...
                    doMerge(ppmCheck == 1) = true;
                else
                    % No point continuing to next variables down the line
                    flag = false;
                end
                
            elseif sum(ppmCheck) == 2
                % See if all three images can be combined at first...
                img1 = double(x(:,:,i) > 0);
                img2 = double(x(:,:,neigh(1)) > 0);
                img3 = double(x(:,:,neigh(2)) > 0);
                
                % Sum together and see if the maximum value is 1
                imgX = img1 + img2 + img3;
                multiCount = sum(imgX(imgX > 1));
                multiPerc = 100 * multiCount / numel(imgX);
                
                if multiPerc < percThreshold
                    % Then images are complementary
                    doMerge = [true true];
                else
                    % Need to determine which is the best of the two...
                    imgY = img1 + img2;
                    imgZ = img1 + img3;
                    
                    maxY = max(imgY(:)) == 1;
                    maxZ = max(imgZ(:)) == 1;
                    
                    sumY = sum(imgY(imgY == 1));
                    sumZ = sum(imgZ(imgZ == 1));
                    
                    mulY = sum(imgY(imgY > 1));
                    mulZ = sum(imgZ(imgZ > 1));
                    
                    perY = 100 * mulY / numel(imgY) < percThreshold;
                    perZ = 100 * mulZ / numel(imgZ) < percThreshold;
                    
                    if sumY >= sumZ && maxY
                        % Then use neigh(1)
                        doMerge(1) = true;                        
                    elseif sumZ > sumY && maxZ
                        % Then use neigh(2)
                        doMerge(2) = true;
                        
                    elseif sumY >= sumZ && perY
                        doMerge(1) = true;                        
                    elseif sumZ > sumY && perZ
                        doMerge(2) = true;
                    else
                        % Turns out that neither were complementary...
                        % Set flag to false to halt the search
                        flag = false;
                    end
                end    
                
            else
                % Don't know how this could happen
                error('Another unexpected exception');
            end
            
            % Now would be the time to deal with the merged variables...
            %doMerge
            %neigh
            %i
            
            if sum(doMerge) > 0
                picked = neigh(doMerge == 1);
                
                % Combine all into the suitable place...
                x(:,:,i) = nansum(x(:,:,[i picked]),3);
                x(:,:,picked) = 0;
                mz(picked) = NaN;
                
                isMerged(picked',1) = true;
                
            else
                
                % Do what?
                
            end
         
            
        end
        
           
        
    else
        % Do nothing?
    end
    
    
end
        
  
fx = isnan(mz);

m2 = mz(~fx);
x2 = x(:,:,~fx);

        
        


return


% Create a figure to show off the results
fig.fig = figure;
fig.ax(1) = subplot(1,3,1); 
fig.im(1) = imagesc(rand(size(x,1),size(x,2)));
fig.ax(2) = subplot(1,3,2); 
fig.im(2) = imagesc(rand(size(x,1),size(x,2)));
fig.ax(3) = subplot(1,3,3); 
fig.im(3) = imagesc(rand(size(x,1),size(x,2)));
linkaxes(fig.ax,'xy');

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

        %set(fig.im(1),'CData',x(:,:,n) > 0);
        %set(fig.im(2),'CData',x(:,:,neigh) > 0);
        %set(fig.im(3),'CData',imgs);
        %drawnow;
        
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


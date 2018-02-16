function [ mz,x ] = combCompNeighbours3(mz,x,ppmTol,flag)
%combCompNeighbours2 - this is a new version which is supposed to be better
%than the first version.  It does this by merging the most likely ones
%before others; i.e. work on highly populated m/z bins first, and ignore
%ones which look most like noise.

percThreshold = 2.5;
verbose = true;

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


% Create a figure to show off the results
if verbose
    fig.fig = figure;
    fig.ax(1) = subplot(2,3,1); 
    fig.im(1) = imagesc(rand(size(x,1),size(x,2)));
    fig.ax(2) = subplot(2,3,2); 
    fig.im(2) = imagesc(rand(size(x,1),size(x,2)));
    fig.ax(3) = subplot(2,3,3); 
    fig.im(3) = imagesc(rand(size(x,1),size(x,2)));
        
    fig.ax(4) = subplot(2,3,4); 
    fig.im(4) = imagesc(rand(size(x,1),size(x,2)));
    fig.ax(5) = subplot(2,3,5); 
    fig.im(5) = imagesc(rand(size(x,1),size(x,2)));
    fig.ax(6) = subplot(2,3,6); 
    fig.im(6) = imagesc(rand(size(x,1),size(x,2)));

    linkaxes(fig.ax,'xy');
end

% Now with that information we can rank the  variables which are better and
% which should be merged first... Let's sort on mean intensity...
[srt,idx] = sortrows(info,-3);

% So could just work down the list of 'srt' to determine which variables
% are potentially merge-able with each of the most populous.  Start with
% one variable either side (within ppm tolerance) and extend only if
% complementary.  If n-1 and n+1 don't work together, then pick the one
% with the best fit.  Don't forget to set merged variables to zero to
% prevent them being used twice.
isMerged = false(numV,1);

% Need something to save all of the statistics to review at the end
all = zeros(numV,6);

for n = 1:numV
    
    % Always refere to its true index, rather than the sorted one
    i = srt(n,1);
    
    % Make sure that this variable hasn't previously been merged
    if isMerged(i,1)
        continue;
    end
    
    % Holes / common pixels
    
    %[~] = holes(x(:,:,i),x(:,:,i+1));
    
    % Otherwise we want to see how well its first two immediate neighbours
    % compare in terms of complementarity...
    values = ones(3,3)*100;
    
    % Center + Left
    if i-1 > 0
        if ~isMerged(i-1,1)
            tmpA = sum(x(:,:,[i i-1]) > 0,3);        
            values(1,1) = overlap(tmpA);

            % Determine holes in the image
            [values(1,2),values(1,3)] = holes(x(:,:,i),x(:,:,i-1));
        end        
    end
    
    % Center + Right
    if i+1 <= numV 
        if ~isMerged(i+1)
            tmpB = sum(x(:,:,[i i+1]) > 0,3);
            values(2,1) = overlap(tmpB);
        
            [values(2,2),values(2,3)] = holes(x(:,:,i),x(:,:,i+1));
        end
    end
    
    % Center + Left + Right
    if i-1 > 0 && i+1 <= numV 
        if ~isMerged(i-1) && ~isMerged(i+1)
            tmpC = sum(x(:,:,i-1:i+1) > 0,3);
            values(3,1) = overlap(tmpC);
        
            [values(3,2),values(3,3)] = holes(x(:,:,i),sum(x(:,:,[i-1 i+1]),3));
        end
    end
    
    % What are the intensities of the points that overlap? Is there a large
    % differnce?
    
    
    %values
    
    if verbose
        set(fig.im(1),'CData',tmpA); %title(fig.im(1),'Center + Left');
        set(fig.im(2),'CData',tmpC); %title(fig.im(2),'Center + Both');
        set(fig.im(3),'CData',tmpB); %title(fig.im(3),'Center + Right');
        
        set(fig.im(4),'CData',x(:,:,i-1) > 0); %title(fig.im(1),'Center + Left');
        set(fig.im(5),'CData',x(:,:,i)   > 0); %title(fig.im(2),'Center + Both');
        set(fig.im(6),'CData',x(:,:,i+1) > 0); %title(fig.im(3),'Center + Right');

        
        drawnow;
    end
    
    
    % Save the relevant statistical results
    all(i,:) = [values(:,2)' values(1:2,3)' values(3,1)];
    
    % So now here we need to decide if the images are worth merging.  We
    % should be looking for a low pixel overlap percentage and a large
    % increase in the number of non-overlapping pixels.
    %values
    
    if values(3,3) < percThreshold && values(3,2) > 20
        % Then merge both images together as they have low overlaps and
        % provide a good addition to the image
        isMerged(i-1:i+1,1) = true;
        
        % Determine weighted average for m/z values
        [mz(i)] = wghtMZ(mz(i-1:i+1),info(i-1:i+1,3));
        mz([i-1 i+1]) = NaN;
        
        % Take the maximum values over the image
        x(:,:,i) = max(x(:,:,[i-1:i+1]),[],3);
        x(:,:,[i-1 i+1]) = 0;
    
    elseif values(1,2) > values(2,2) && values(1,3) < percThreshold
        % Then merge the left hand image only
        isMerged(i-1:i,1) = true;
        
        % Determine weighted average for m/z values
        [mz(i)] = wghtMZ(mz(i-1:i),info(i-1:i,3));
        mz([i-1]) = NaN;
        
        % Take the maximum values over the image
        x(:,:,i) = max(x(:,:,[i-1:i]),[],3);
        x(:,:,i-1) = 0;

        
    elseif values(2,2) > values(1,2) && values(2,3) < percThreshold
        % Then merge the right hand image only
        isMerged(i:i+1,1) = true;
        
        % Determine weighted average for m/z values
        [mz(i)] = wghtMZ(mz(i:i+1),info(i:i+1,3));
        mz([i+1]) = NaN;
        
        % Take the maximum values over the image
        x(:,:,i) = max(x(:,:,[i:i+1]),[],3);
        x(:,:,i+1) = 0;
        
    else
        % Do no merging here...
        
    end
        
    
    
end

% Trim out the dross
fx = isnan(mz);
mz = mz(~fx);
x = x(:,:,~fx);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ol] = overlap(img)
% Determine overlap percentage

f2 = img(:) == 2;

ol = 100 * sum(f2) / numel(img);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hfs,cps] = holes(a,b)
% Determine how many holes in a that b fills
% Determine how many pixels that both a and b have

%f1 = a > 0;
%f2 = b > 0;

% Holes filled
hf = a == 0 & b > 0;

% Common pixels
cp = a > 0 & b > 0;

hfs = 100 * sum(hf(:)) / numel(a);
cps = 100 * sum(cp(:)) / numel(a);

% For the overlapped/common pixels, we want to see what their intensities
% are like
%ia = a(cp);
%ib = b(cp);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = wghtMZ(mz,mv)
% Determine new mz based on weighted average

wtm = sum(mz(:) .* mv(:)) / sum(mv(:));
new = wtm;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
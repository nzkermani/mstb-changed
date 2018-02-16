function [ data,origMZ ] = rawH5extract(file,pks)
%rawH5extract - using the peak list, extract the images from a rawH5 file

tic

intMethod = 'max';

% How many peaks are there?
numPks = numel(pks.mz);

% How big is the file?
totInt = h5read(file,'/totPts');
numPix = numel(totInt);

% Image xy coordinates
xy = h5read(file,'/xy2D');
xy = reshape(xy,[size(totInt,1)*size(totInt,2) 1]);

% Create a matrix for storage. Let's do it in two 2D to make life easier to
% start with. Perhaps also possible in sparse format
data = zeros(numPix,numPks);
origMZ = zeros(numPix,numPks);

switch intMethod
    case 'max'
        fInt = @(x) max(x);
    case 'sum'
        fInt = @(x) sum(x);
end

% Create a waitbar
wb = waitbar(0,file);

% Loop through each pixel...
for n = 1:numPix
        
    % Determine the xy index so that we can reshape the image
    i = xy(n);

    % Scan name
    sn = ['/raw/scan/' int2str(i)];
        
    % Read spectrum...
    tmp = h5read(file,sn);
    
    % Loop through each of the peaks
    for r = 1:numPks
        
        % What is the peak span?
        span = pks.span(r,:) - pks.mz(r);
        wdth = min(abs(span));
        
        % What are the crucial dimensions of the peak?
        locn = [pks.mz(r) - wdth pks.mz(r) + wdth];
        
        % Find either the largest intensity within this boundary, or sum
        % all points within it...
        fx = tmp(:,1) >= locn(1) & tmp(:,1) <= locn(2);
        
        % Or should we be identifying the nearest peak instead?
        
        if sum(fx) > 0
            data(n,r) = fInt(tmp(fx,2));
            
            % Determine the raw m/z value
            ppp = tmp(fx,:);
            [~,tmpIdx] = max(ppp(:,2));
            origMZ(n,r) = ppp(tmpIdx,1);
        end
        
    end
    
    waitbar(n/numPix,wb);
end

toc

data = reshape(data,[size(totInt,1) size(totInt,2) numPks]);
origMZ = reshape(origMZ,[size(totInt,1) size(totInt,2) numPks]);

delete(wb);

end


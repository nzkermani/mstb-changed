function [ recal ] = globalRecal(data,flag)
%globalRecal - perform global recalibration of an imzML file using a few
%sentinel ions

if nargin == 1
    flag = false;
end

% Define the ions for assessment
ions = [255.2329 281.2485 465.3043 536.5047 744.5548 885.5498]';
hw = 15;

% Extract images
[mz,~] = getImages(data,ions,hw);

% So now that we have the images, we need to determine the median m/z value
% for each ion
[meds] = getMedianMZ(mz);

% Determine gross ppm deviation
dev = 1e6 * (meds-ions) ./ ions;

if flag
    recal = dev;
    return
end

% Determine a linear fit for this data - using ppm scales should help this
% considerably. In the case of totally abnormal values, then we need to
% flag this, and consider ditching the middle ones if they seem more
% outlying than the others
[pp] = fitDeviations(ions,dev);

% So now we need to perform recalibration on each scan...
[recal] = correctScanByScan(data,pp);

return

% Go back and extract the images for the recalibrated data, and then
% determine the ppm shift again to check it looks good
[rm,rs] = getImages(recal,ions,3);
[rmeds] = getMedianMZ(rm);
dev = 1e6 * (rmeds-ions) ./ ions;

% Make a nice image...
diagnosticImage(rm,rs,ions);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mz,im] = getImages(data,ions,hw)

% Image storage
sz = size(data);
mz = NaN(sz(1),sz(2),numel(ions));
im = NaN(sz(1),sz(2),numel(ions));

% Generate ion images for these ions, taking the largest peak within the
% window. Then we will use the median values to determine a calibration
% function
for p = 1:size(data,1)
    for q = 1:size(data,2)
        
        % Loop through each ion
        for n = 1:numel(ions)
            
            % Save spectrum
            tmp = data{p,q};
            if isempty(tmp)
                continue;
            end                
    
            % Find largest signal
            df = hw * ions(n) / 1e6;
            fx = tmp(:,1) > (ions(n) - df) & tmp(:,1) < (ions(n) + df);
            
            % Add to image
            if sum(fx) > 0
                fx = fx ./ tmp(:,2);                
                [~,b] = max(fx);    
                
                mz(p,q,n) = tmp(b,1);
                im(p,q,n) = tmp(b,2);           
            end
        end
    end
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [meds] = getMedianMZ(mz)
% Determine the median m/z values for each of the ions

numI = size(mz,3);
meds = NaN(numI,1);

for n = 1:numI
    
    tmp = mz(:,:,n);
    meds(n,1) = nanmedian(tmp(:));
    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = fitDeviations(ions,dev)
% Linear fit

p = polyfit(ions,dev,1);

return

% Extrapolate over normal m/z range for figure vis purposes
mz = 0:1:1000;
proj = polyval(p,mz);
figure; hold on;

plot(mz,proj);
stem(ions,dev,'r');


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = correctScanByScan(data,pp)
% Recalibrate each scan sequentially

% Now is time to recalibrate all of the scans...
for p = 1:size(data,1)
   
    % Each of the row's pixels
    for q = 1:size(data,2)
        
        % Spectrum for scan p,q
        tmp = data{p,q};
        if isempty(tmp)
            continue;
        end
        
        % Apply function to determine ppm shift
        shift = polyval(pp,tmp(:,1));
        
        % Determine the offset in m/z units (i.e. convert from ppm)
        offset = tmp(:,1) .* shift / 1e6;
        
        % New mz values...
        newMZ = tmp(:,1) - offset;
        
        % Save...
        data{p,q}(:,1) = newMZ;
        
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diagnosticImage(rm,rs,ions)
% Figure to show ion images of recalibrated m/z and intensities

numI = numel(ions);

figure('Units','normalized','Position',[0 0.5 1 0.5]);
ax = zeros(2,numI);

for n = 1:numI
    
    ax(1,n) = subplot(2,numI,n);
    imagesc(rs(:,:,n));
    colorbar;
    
    ax(2,n) = subplot(2,numI,n+numI);
    tmp = 1e6 * (rm(:,:,n)-ions(n)) / ions(n);
    imagesc(tmp);
    colorbar;
    caxis([-3 3]);

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


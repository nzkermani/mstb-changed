function [ recal ] = imageRecal(data)
%imageRecal - recalibration of imzml data

% Define the ions for assessment
ions = [255.2329 281.2485 465.3043 536.5047 744.5548 885.5498];
hw = 20;

% Create a blank matrix in which to store the mz values
numI = numel(ions);
vals = NaN([size(data) numI]);

for p = 1:size(data,1)    
    for q = 1:size(data,2)
        if ~isempty(data{p,q})
            vals(p,q,:) = scanRecal(data{p,q}(:,1),data{p,q}(:,2),ions,hw);
        end
    end
end

% Check ppm deviations...
ppms = bsxfun(@minus,vals,reshape(ions,[1 1 numel(ions)]));
ppms = 1e6 * bsxfun(@rdivide,ppms,reshape(ions,[1 1 numel(ions)]));

% Somewhere to store the polynomial coefficients
ord = 1;
coeff = NaN(size(data,1),ord+1);

% Can we determine the mean m/z values across each row?
for p = 1:size(data,1)
    
    % Determine the median 
    md = squeeze(nanmedian(ppms(p,:,:)));
    
    % Only can do it if we have 3 points? Pref first and last?
    numP = sum(~isnan(md));
    endP = ~isnan(md(1)) & ~isnan(md(end));
    
    % Do?
    if numP >= 3 && endP
        
        fx = ~isnan(md);
        
        % Fit a 1st order polynomial
        coeff(p,:) = polyfit(ions(fx)',md(fx),1);

        % Expand over the full m/z range
        %y1 = polyval(,100:0.1:1000);        
        
    else
        % Do nothing else...
        continue;
    end
    
end

% Now that we have all of the coefficients, we can guess the coefficients
% for the empty pixels by using the average value over the neighbouring ±3
% scans
wn = 3;
fx = find(isnan(coeff(:,1)));
newC = coeff;
for n = 1:numel(fx)
    
    range = [max([1 fx(n)-wn]) min([fx(n)+wn size(coeff,1)])];
    
    tmp = nanmean(coeff(range,:),1);
    
    % If still no value use the global mean values
    if isnan(tmp(1))
        tmp = nanmean(coeff,1);
    end
    
    % Save to new variable
    newC(fx(n),:) = tmp;
    
end

% Now is time to recalibrate all of the scans...
recal = cell(size(data));
for p = 1:size(data,1)
    
    % This is the polynomial function to be used
    pp = newC(n,:);
    
    % Each of the row's pixels
    for q = 1:size(data,2)
        
        oldMZ = data{p,q}(:,1);
        
        % Apply function
        ppmShift = polyval(pp,oldMZ);        
        offset = oldMZ .* ppmShift / 1e6;        
        newMZ = oldMZ + offset;
        
        % Save...
        recal{p,q} = [newMZ data{p,q}(:,2)];
        
    end
    
end
    
    


end


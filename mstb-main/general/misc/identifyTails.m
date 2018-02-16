function [final] = identifyTails(x,y,g,flag)
%identifyTails - find variables at the base of large peaks that are
%artifacts. Try to remove them from the spectra

% For each group (g) find the average spectrum.
% For large peaks determine if it has peaks within a window of m/z|int

% For plotting
if nargin == 3
    flag = false;
end

% Max peaks to look at?
maxP = 10; % max number of peaks to look at
mhw = 0.20; % m/z half width
prc = 0.01; % intensity fraction

% Group information
[unqG,~,indG] = unique(g);
numG = numel(unqG);
numV = numel(x);

% Create a vector to see which variables are to be removed
ditch = zeros(numV,numG);
means = zeros(numV,numG);

for n = 1:numG
    
    fx = indG == n;
    
    % Mean spectrum
    mn = nanmean(y(fx,:),1);
    means(:,n) = mn';
    
    % Find the maxP most intense peaks
    [srt,idx] = sort(mn,'descend');
    idx = idx(1:maxP);
    
    for r = 1:maxP
        i = idx(r);
        
        % These are the variables within the m/z window with intensities
        % less than the fraction of the super intense peak
        xfail = x > x(i)-mhw & x < x(i)+mhw & mn < mn(i)*prc;
        
        % To ensure that these are proper tails, we have to ensure that
        % there are about 3 tailing peaks either side of the main peak.
        chkL = sum(xfail(1:i));
        chkR = sum(xfail(i:end));        
        if chkL >= 3 && chkR >= 3       
            ditch(:,n) = ditch(:,n) + double(xfail');
        else
            % Do nothing, as there aren't enough tails for this to be a
            % likely set of tailing peaks
            %disp('not tails');
        end        
        
        %figure; hold on; stem(x,mn,'k'); stem(x(xfail),mn(xfail),'r');
        %xlim([x(i)-(2*mhw) x(i)+(2*mhw)]);
        
    end
    
end

final = sum(ditch,2) > 0;

if ~flag
    return
end

figure; hold on;
stem(x,means);
stem(x(final),-means(final,:),'bs');

% Create the new matrix
x2 = x(~final);
y2 = y(:,~final);
z2 = y(:,final);

% PCA to see if the data is changed a lot
[l1,s1,e1] = princomp(bsxfun(@rdivide,y,nansum(y,2)));
[l2,s2,e2] = princomp(bsxfun(@rdivide,y2,nansum(y2,2)));
[l3,s3,e3] = princomp(bsxfun(@rdivide,z2,nansum(z2,2)));

figure; hold on;

subplot(1,3,1);
scatter(s1(:,1),s1(:,2),80,indG,'o','filled');
title('Original');

subplot(1,3,2);
scatter(s2(:,1),s2(:,2),80,indG,'o','filled');
title('Tails Removed');

subplot(1,3,3);
scatter(s3(:,1),s3(:,2),80,indG,'o','filled');
title('Just the Tails');


end


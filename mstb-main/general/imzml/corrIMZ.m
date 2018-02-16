function [ newMZ,newImg,opts ] = corrIMZ(unqMZ,img,opts)
%corrIMZ - correlate neighbouring ion images

% Options for defaults
if nargin == 2
    opts.occupancy  = 0.005;
    opts.corrThresh = 0.200;
    opts.connect    = 2;        % 1 would include everything
end

% First we reshape from 3D to 2D to make life easier
sz = size(img);
orig = img;
img = reshape(img,[sz(1)*sz(2) sz(3)]);

% Blank space for correlation coefficient
crr = zeros(sz(3),3);

if opts.occupancy < 1
    multFac = sz(1) * sz(2);
else
    multFac = 1;
end

% Correlate each image to its two neighbours, saving the maximum value
notice = round(linspace(1,sz(3),20));
disp('>>');
for n = 1:sz(3)-1
    
    % We need to filter out the zero-zero values (or NaN-NaN values)
    tmp = img(:,[n n+1]);
    
    % This works for the non zero-zero pairs
    f1 = sum(tmp > 0,2) == 2;
    
    % What if there are no points?
    if sum(f1) < 2
        xx = NaN;        
    else    
        % Calculate
        xx = corr(img(f1,n),img(f1,n+1));
    end
    
    % Place in the matrix
    crr(n,1) = max([crr(n,1) xx]);
    crr(n+1,1) = xx;
        
    % How many pairs of data is this made up from?
    crr(n,2) = sum(f1) / multFac;
    
    % Perform the connectivity function to ensure that we have at least two
    % neighbouring pixels...
    conn = bwconncomp(orig(:,:,n) > 0);
    szes = cellfun(@size,conn.PixelIdxList','UniformOutput',false);
    if ~isempty(szes)
        szes = vertcat(szes{:});
        crr(n,3) = max(szes(:,1));
    end
    
    % Percentage complete
    if any(notice == n)
        disp([char(8) int2str(find(notice == n)) ' ']);
    end        
end
    
% I think we can easily exclude peaks that have fewer than 10 points and a
% negative correlation. The rest is unknown, but at the discretion of the
% user perhaps...
ic = crr(:,1) > opts.corrThresh & ...   % above correlation minimum
    crr(:,2) > opts.occupancy & ...     % at least x% pixel occupation
    crr(:,3) >= opts.connect;           % at least y neighbouring pixels

% Expect to find non-isolated peaks
il = [ic(2:end); false];
ih = [false; ic(1:end-1)];

inc = sum([il ic ih],2) > 2 & ic == 1;

inc(find(inc)-1) = true;
inc(find(inc)+1) = true;

% Plot a graph showing everything
% figure; hold on;
% plot(unqMZ(1:sz(3)),crr(:,1),'LineWidth',0.5);
% scatter(unqMZ(inc),crr(inc,1),80,'red','o','filled');

disp(['Keep:' int2str(sum(inc)) ' / ' int2str(sz(3))]);

% An average mass spectrum on which to project the kept variables
avg = nanmean(img,1);

plotXYC(unqMZ,avg,crr(:,1)); hold on;
xlabel(gca,'m/z  ','FontWeight','bold','FontSize',20,'FontAngle','italic');
ylabel(gca,'Intensity (mean)  ','FontWeight','bold','FontSize',20);
cb = colorbar;
ylabel(cb,'Correlation  ','FontWeight','bold','FontSize',20);
set(gca,'FontSize',16);
box on;
scatter(unqMZ(inc),avg(inc),80,'red','o','filled');

newImg = img(:,inc);
newImg = reshape(newImg,[sz(1) sz(2) sum(inc)]);
newMZ  = unqMZ(inc);

end

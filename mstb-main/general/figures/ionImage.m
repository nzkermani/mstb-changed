function [ output_args ] = ionImage(x,r,g,b)
%ionImage - using the matrix of x, make an ion image using the supplied RGB
%channels...

sz = size(x);

img = zeros(sz(1),sz(2),3);

for n = 1:3
    
    % Define the indices
    if n == 1
        i = r;
    elseif n == 2
        i = g;
    else
        i = b;
    end
    
    % Just blank
    if isempty(i)
        continue;
    end
    
    % Sum all ions together
    tmp = nansum(x(:,:,i),3);
    
    % Set super high intensities a little lower
    p75 = prctile(tmp(tmp > 0),95);
    tmp(tmp > p75) = p75;
    
    % Smooth?
    filt = fspecial('gaussian',5);
    tmp = filter2(filt,tmp);
   
    
    img(:,:,n) = tmp;
    
end

% Scale parts betweeen 0 and 1
img = imScale(img);

% Plot the image
fig.fig = figure('Name','Ion Image',...
    'Units','normalized',...
    'Position',[0.5418 0.2146 0.3879 0.6903]);

hold on;

% Draw the axes...
fig.ax1 = axes('Parent',fig.fig,...
    'Units','normalized',...
    'Position',[0 0 1 1],...
    'XTick',[],...
    'YTick',[]);

imagesc(img,'Parent',fig.ax1);


%'Units','normalized','Position',[0 0 1 1],...
    


end


function interpExp(data)
%interpExp - exploration of the benefits / problems of interpolation of
%profile mode MS data
%
%
% Run this using Jocelyn's OrbData1:
% /Volumes/Data/Data/Misc/Jocelyn Sprayer Data/Batch1_OrbData.mat
%
% James McKenzie, 2015.

% Define a couple of custom interpolation ranges
ppmRes = [2.5 10];      % percieved to be a good couple of values
dalRes = [0.001 0.01];  % the classic ranges
mzR    = [200 1500];    % m/z range over which to care

% Now (actually however fill) we will determine the polynomial fit for the
% dataset in question. NB this function won't work in its current form for
% other data formats.
[fit] = detPPMdiffs(data,true)
%fit = [-0.000001375942518 0.005312507060135 1.399428855228931];
fit = [0.000000000828250 -0.000003169555291 0.006424315201589 1.221522478597406];

fit(end) = fit(end) * 2

% What kind of interpolation method to perform? (linear/pchip)
intMet = {'linear'};

% Now create the m/z vectors
[vec] = mzVectors(mzR,ppmRes,dalRes,fit);
numV = size(vec,2);

% How many observations in the dataset?
numO = size(data,2);

% So can run through the entire dataset, performing different forms of
% interpolation and stuff...


for n = 1%1:numO
    
    % How many scans in the observation?
    numS = 5;%size(data{1},2);
    
    % Can we interpolate this data?    
    resSt(numS,numV) = struct('mz',[],'sp',[],'lab',[]);
    
    for r = 1:numS
                
        % Define the mz mask, generically useful 
        msk = data{n}{r}(:,1) > mzR(1) & data{n}{r}(:,1) < mzR(2);
        
        
        % Now do the interpolation for each of the resolutions/vectors as
        % calculated above        
        for v = 1:numV
            
            % Instead of that stuff, let's shift the variables a little 
            % each time to see how the interpolation handles little shifts
            %interpShift(data{n}{r}(msk,1),data{n}{r}(msk,2),vec(v).mz,vec(v).label);

        
            % Call the interpolation
            resSt(r,v).sp = doInterp(vec(v).mz,...    % new m/z vector
                    data{n}{r}(msk,1),...           % old m/z vector
                    data{n}{r}(msk,2),...           % old intensities
                    intMet{1});                     % interp method
            
            % Do we need to save the mz vector again?
            resSt(r,v).mz = vec(v).mz;
                        
            % Create a label for it
            resSt(r,v).lab = vec(v).label;
        
        end
        
        % Here once we've done the calculations we can plot the data in a
        % good way that shows the scumsum, interpolation, and high and low
        % m/z peaks for interp1 comparison
        figComp(resSt(r,:),data{n}{r}(msk,1),data{n}{r}(msk,2));
                       
        % Here is where we instigate the accumarray binning method. This is
        % notionally easier for fixed width data as it is periodic
        % [aa,bb] = doAccArr(data{n}{r}(msk,1),data{n}{r}(msk,2),0.001,@sum);
        % [cc,dd] = doAccArr(data{n}{r}(msk,1),data{n}{r}(msk,2),0.001,@mean);
        % 
        % figure; hold on;
        % plot(data{n}{r}(msk,1),data{n}{r}(msk,2),'-ok');
        % plot(aa,bb,'-ob');
        % plot(cc,dd,'-or');
               

    end
    
    % Determine the variances over the r scans    
    [allv] = detVar(resSt,vec);
    

    
    
    
end










end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = doInterp(interpMZ,mz,sp,method)
% Perform interpolation according to the specified method

y = interp1(mz,sp,interpMZ,method);

y(isnan(y)) = 0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec] = mzVectors(mzR,ppmRes,dalRes,fit)
% Create the various vectors for the analyses

vec = struct('mz',[],'res',[],'type',[],'label',[]);

allRes = [ppmRes dalRes];
for n = 1:numel(allRes)
    
    if allRes(n) >= 1
        vec(n).mz = ppmVector(mzR(1),mzR(2),allRes(n))';                
        vec(n).type = 'ppm';
        vec(n).label = {[num2str(allRes(n)) ' ppm']};            
    else
        vec(n).mz = mzR(1):allRes(n):mzR(2);
        vec(n).type = 'dal';
        vec(n).label = {[num2str(allRes(n)) ' Da']};    
    end    
    vec(n).res = allRes(n);
end

% If we have provided a 'fit' polynomial then do that one
if nargin == 4
    vec(n+1).mz = ppmVector(mzR(1),mzR(2),fit)';
    vec(n+1).type = {'super'};
    vec(n+1).label = {'variable ppm'};
    vec(n+1).res = fit;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new,ac] = doAccArr(mz,sp,res,flag)
% Perform accumulated binning...

% First we round the mz values
mzN = round(mz / res);

minN = mzN(1);
mzN = mzN - mzN(1) + 1;

% Accumulate...
ac = accumarray(mzN,sp,[],flag);
    
% Generate the new MZ scale
new = mzN(1):1:mzN(end);
new = new + minN - 1;
new = new * res;

fx = ac > 0;

ac = ac(fx);
new = new(fx);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = scumsum(x)
% Scale the cumsum to 100

y = cumsum(x);
y = 100 * y / y(end);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figComp(res,origX,origY)
% Draw a figure that beautifully expresses the interpolated data and how
% and where and why it fails...

fig.fig = figure('Units','normalized','Position',[0 0 1 1]);

% Draw some subplots
%fig.ax(2) = subplot(2,2,2);
%fig.ax(4) = subplot(2,2,4);

% In the first one we do the scumsums
fig.ax(1) = subplot(2,2,[1 2]); hold on;

plot(origX,scumsum(origY),...
    'Color',[0 0 0],...
    'LineWidth',3);

[~,numV] = size(res);
cols = jet(numV);
if numV == 5
    cols(4,:) = [225 21 232] / 256;
elseif numV == 4
    cols(3,:) = [225 21 232] / 256;
end

hand = zeros(numV,2);
for n = 1:numV    
    % Plot!
    hand(n,1) = plot(res(n).mz,scumsum(res(n).sp),...
        'Color',cols(n,:),...
        'LineWidth',2);    
end

% Legend?
legend(hand(:,1),vertcat(res.lab),'Location','NorthWest');

% Labels
ylabel('Scaled cumulative sum','FontSize',16,'FontWeight','bold');
xlabel('m/z','FontAngle','italic','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14);


% How about two plots showing low/high m/z ranges...
%fig.ax(3) = subplot(2,2,3); hold on;
mzrange = [269.081 269.094; 885.507 885.596];

for r = 1:size(mzrange,2)
    
    % Define the axes
    if r == 1
        fig.ax(3) = subplot(2,2,3); hold on;
    else
        fig.ax(4) = subplot(2,2,4); hold on;
    end
    
    sulab = cell(numV+1,1);
    
    % Original data...
    msk = origX > mzrange(r,1) & origX < mzrange(r,2);
    tot = stem(origX(msk),origY(msk),...
        'Color',[0 0 0],...
        'LineWidth',2);
    sulab{numV+1,1} = [int2str(sum(msk)) ' points'];
        
    mzzz = origX(msk);
    labb = ['m/z range = ' sprintf('%0.4f',max(mzzz)-min(mzzz))];
    xlabel(labb,'FontSize',16,'FontWeight','bold');
    ylabel('Intensity','FontSize',16,'FontWeight','bold');
    
    % Now consider looping through the data...
    
    for n = 1:numV
        % Plot!
        msk = res(n).mz > mzrange(r,1) & res(n).mz < mzrange(r,2);
        msk2= res(n).sp > 0;
        msk = (msk .* msk2) == 1;
        
        sulab{n,1} = [int2str(sum(msk)) ' points'];
        
        if sum(msk) == 1
            ff = find(msk == 1);
            msk(ff-1) = true;
            msk(ff+1) = true;
        end
        
        try
            hand(n,2) = plot(res(n).mz(msk),res(n).sp(msk),'--',...
                'Color',cols(n,:),...
                'LineWidth',4,...
                'MarkerFaceColor',cols(n,:),...
                'MarkerEdgeColor',cols(n,:)); 
        end
    end
    
    legend([hand(:,2)' tot],sulab);
    xlim(mzrange(r,:));
    set(gca,'FontSize',14);

end






end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Y,P,shftOpt] = msalignJSM(X,Y,P,varargin)
%MSALIGN signal calibration and alignment by reference peaks
%
%   YOUT = MSALIGN(X,Y,P) Aligns a signal with peaks (Y) to a set of
%   reference peaks (P) by scaling and shifting the domain (X) such that
%   the cross-correlation between the input signal and a synthetic target
%   signal is maximum. The synthetic target signal is built with Gaussian
%   pulses centered at the locations specified by the vector P. After
%   obtaining the new X scale, MSALIGN calculates YOUT by shape-preserving 
%   piecewise cubic interpolation of the shifted input signal to the
%   original X vector. 
%
%   X and Y are column vectors where paired values represent points in the
%   signal. Y can be a matrix with several signals, all sharing the same X
%   scale. Units in the X scale (separation units or s.u.) may quantify
%   wavelength, frequency, distance, time or m/z depending on the type of
%   instrument that generates the signal.
%
%   Notes:
%   1) MSALIGN uses an iterative grid search until it finds the best scale
%   and shift factors for every signal. 
%   2) Better calibration is achieved when more reference peaks are known.
%   3) When only one reference peak is given the algorithm only aligns the
%   signal with constant shifts and does not rescale. 
%   4) MSALIGN also allows to achieve multiple alignment of signals by
%   providing only a coarse estimation of a set of common peaks that
%   appears in every signal, see the 'GROUP' option below.  
%
%   MSALIGN(...,'RESCALING',FALSE) turns off the rescaling of X, the output
%   signal is only aligned to the reference peaks by using constant shifts.
%   By default a rescaling factor is estimated by MSALIGN unless one
%   reference peak (P) is used.
%
%   MSALIGN(...,'WEIGHTS',W) sets the weighting for every peak in P. W is a
%   vector of the same size as P and it contains relative weights with
%   non-negative values. The default is W = ones(size(P)). Increase the
%   relative weightings to emphasize small reference peaks in the signal.
%
%   MSALIGN(...,'MAXSHIFT',MS) sets the lower and upper limits for allowable
%   shifts in range for any of the peaks, the default is [-100 100] s.u.
%   Use these values to tune the robustness of the algorithm. Ideally you
%   should only try to correct small shifts; therefore, these bounds should
%   be small. However, by increasing the limits is possible to correct
%   larger shifts at the expense of picking incorrect peaks to be aligned.
%
%   [YOUT, POUT] = MSALIGN(...,'GROUP',TRUE) updates the location of the
%   reference peaks (P) such that the actual movement of the peaks in the
%   signals is minimum. POUT contains the updated peak locations. The
%   default is FALSE.
%
%   MSALIGN(...,'SHOWPLOT',SP) plots the original and the aligned signals
%   over the P markers. When SP is TRUE, the first signal in Y is used. If
%   MSALIGN is called without output arguments, a plot will be shown unless
%   SP is FALSE. SP can also contain an index to one of the signals in Y. 
%
%   Advanced algorithmic parameters:
%
%   MSALIGN(...,'WIDTHOFPULSES',WP) sets the width WP (in s.u.) of the
%   Gaussian pulses used to build the correlating synthetic signal. WP is
%   the point where the Gaussian pulse reaches 60.65% of its maximum. The
%   default is 10, which means a constant width for all Gaussian pulses is
%   set at 10. WP can also be a function handle. The referenced function is
%   evaluated at each X value to compute a variable width for the pulses.
%   Its evaluation should give reasonable values (i.e., WP < max(abs(R) and
%   WP >0); otherwise, MSALIGN errors out.
%
%   Note: Tuning the spread of the Gaussian pulses controls a tradeoff
%   between robustness (wider pulses) and precision (narrower pulses);
%   however, it is unrelated to the shape of the observed peaks in the
%   signal. 
%
%   MSALIGN(...,'WINDOWSIZERATIO',WSR) WSR is a scalar value that
%   determines the size of the windows around every putative peak; the
%   synthetic signal is correlated to the input signal only within these
%   regions, saving computation time. The size of the window is given by
%   WP * WSR (in s.u.). The default is 2.5; i.e., at the limits of the
%   window, the Gaussian pulses have a value of 4.39% of their maximum. 
%
%   MSALIGN(...,'ITERATIONS',I) sets the number of refining iterations. At
%   every iteration, the search grid is scaled down to improve the
%   estimates. The default is 5.
%
%   MSALIGN(...,'GRIDSTEPS',GS) sets the number of steps for the search
%   grid; i.e., at every iteration the search area is divided by GS^2.
%   The default is 20. 
% 
%   MSALIGN(...,'SEARCHSPACE',SS) sets the type of the search space. The
%   default is a 'regular' grid. 'latin' uses a random latin hypercube with
%   GS^2 samples.
%
%  
%   Example: 
%
%       load sample_lo_res
%       % Select some reference peaks.
%       P = [3991.4 4598 7964 9160];
%       msheatmap(MZ_lo_res,Y_lo_res,'markers',P,'range',[3000 10000])
%       title('before alignment')
%
%       % Align all the spectrograms in Y_lo_res to the reference peaks.
%       YA = msalign(MZ_lo_res,Y_lo_res,P);
%
%       msheatmap(MZ_lo_res,YA,'markers',P,'range',[3000 10000])
%       title('after alignment')
%
%       % Repeat the alignment now specifying weights for the reference
%       % peaks.
%       P = [3991.4 4598 7964 9160];
%       W = [60 100 60 100];
%       YA = msalign(MZ_lo_res,Y_lo_res,P,'weights',W);
%
%       msheatmap(MZ_lo_res,YA,'markers',P,'range',[3000 10000])
%       title('after alignment with weights')
%
%   See also MSBACKADJ, MSHEATMAP, MSLOWESS, MSNORM, MSPREPRODEMO,
%   MSRESAMPLE, MSSGOLAY, MSVIEWER. 

%   Copyright 2003-2008 The MathWorks, Inc.


% check inputs
bioinfochecknargin(nargin,3,mfilename);
% set defaults
onlyshift = false;
gaussianWidth = 10;  % std dev of the Gaussian pulses (in X)
gaussianRatio = 2.5; % sets the width of the windows use at every pulse
gaussianResol = 100; % resolution of every Gaussian pulse (number of points) 
iterations = 5;      % increase to improve accuracy
gridSteps = 20;      % size of grid for exhaustive search      
shiftRange = [-100 100];     % initial shift range
searchSpaceType = 'regular'; % type of search space
groupAlign = false;
if nargout == 0
    plotId = 1; 
else
    plotId = 0;
end

% validate required inputs 
if ~isnumeric(Y) || ~isreal(Y) 
   error(message('bioinfo:msalign:IntensityNotNumericAndReal')) 
end
if ~isnumeric(X) || ~isreal(X) || ~isvector(X)
   error(message('bioinfo:msalign:XNotNumericAndReal')) 
end
if size(X,1) ~= size(Y,1)
   error(message('bioinfo:msalign:NotEqualNumberOfSamples'))
end
numSignals = size(Y,2);
if ~isnumeric(P) || ~isreal(P) || ~isvector(P)
   error(message('bioinfo:msalign:PValsNotNumericAndReal')) 
end

P=P(:);
numPeaks = numel(P);
if numPeaks==1
    onlyshift = true;
end
W = ones(size(P));  %default weights


% get input arguments
if  nargin > 3
    if rem(nargin,2) == 0
        error(message('bioinfo:msalign:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'weights','range','widthofpulses',...
              'windowsizeratio','size','iterations','gridsteps',...
              'steps','searchspace','space','group','showplot',...
              'rescaling','maxshift'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:msalign:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:msalign:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % 'weights'
                    W =  pval(:);
                    if numel(W)~=numPeaks
                        error(message('bioinfo:msalign:IncorrectWeights'));
                    end
                case {2,14} % range maxshift
                    shiftRange = pval(:)';
                    if (numel(shiftRange)~=2) || (diff(shiftRange)<=0 )
                        error(message('bioinfo:msalign:IncorrectShiftRange'));
                    end
                case 3 % width
                    gaussianWidth = pval;
                case {4,5} % window size
                    if (~isscalar(pval) || pval<=1)
                        error(message('bioinfo:msalign:IncorrectWindowRatio'))
                    end
                    gaussianRatio = pval;
                case 6 % iterations
                    if (~isscalar(pval) || pval<=0)
                        error(message('bioinfo:msalign:IncorrectPrecision'))
                    end
                    iterations = pval;
                case {7,8} %  grid steps
                    if (~isscalar(pval) || pval<=0)
                        error(message('bioinfo:msalign:IncorrectGrid'))
                    end
                    gridSteps = pval; 
                case {9,10} % search space
                    searchSpaceTypes = {'latinhypercube','regular'};
                    searchSpaceType = strmatch(lower(pval),searchSpaceTypes); 
                    if isempty(searchSpaceType) 
                        error(message('bioinfo:msalign:NotValidSearchSpaceType'))
                    end
                    searchSpaceType = searchSpaceTypes{searchSpaceType};
                case 11 % group
                    groupAlign = bioinfoprivate.opttf(pval,okargs{k},mfilename);
                case 12 % show
                    if bioinfoprivate.opttf(pval) 
                        if isnumeric(pval)
                            if isscalar(pval)
                                plotId = double(pval); 
                            else
                                plotId = 1;
                                warning(message('bioinfo:msalign:SPNoScalar'))
                            end
                        else
                            plotId = 1;
                        end
                    else
                        plotId = 0;
                    end
                case 13 % rescaling
                    onlyshift = ~bioinfoprivate.opttf(pval,okargs{k},mfilename);
            end
        end
    end
end

if (plotId~=0) && ~any((1:numSignals)==plotId)
    warning(message('bioinfo:msalign:InvalidPlotIndex'))
end

% change scalar to function handler
if isnumeric(gaussianWidth)   
    gaussianWidth = @(x) gaussianWidth;   
end

% check that values for GaussianWidth are valid
G = zeros(numPeaks,1);
for i = 1:numPeaks
    G(i) = gaussianWidth(P(i));
    if G(i)<=0 || G(i)>max(abs(shiftRange))
        error(message('bioinfo:msalign:InvalidWidths'))
    end
end

% set the synthetic target signal
corr_sig_x = zeros(gaussianResol+1,numel(P)); 
corr_sig_y = zeros(gaussianResol+1,numel(P)); 
for i = 1:numPeaks
    leftL = P(i)-gaussianRatio*G(i);
    rightL = P(i)+gaussianRatio*G(i);
    corr_sig_x(:,i) = (leftL+(0:gaussianResol)*(rightL-leftL)/gaussianResol)';
    corr_sig_y(:,i) = W(i)*exp(-((corr_sig_x(:,i)-P(i))/G(i)).^2);
end

corr_sig_l = (gaussianResol+1)*numPeaks;
corr_sig_x = corr_sig_x(:);
corr_sig_y = corr_sig_y(:);

% set reduceRangeFactor to take 5 points of the previous ranges or half of
% the previous range if gridSteps<10
reduceRangeFactor = min(0.5,5/gridSteps);

% set scl such that the maximum peak can shift no more than the
% limits imposed by shft when scaling
scaleRange = 1 + shiftRange/max(P); 

if onlyshift
   scaleRange = [1 1];
end

% allocate space for vectors of optima
sclOpt = zeros(numSignals,1); shftOpt = sclOpt;

% number of points in the search space
searchSize = (gridSteps)^2;

% create the meshgrid only once
switch searchSpaceType
    case 'regular'
        [A,B] = meshgrid((0:gridSteps-1)/(gridSteps-1),...
                         (0:gridSteps-1)/(gridSteps-1));
        srchSpc = repmat([A(:),B(:)],1,iterations);
    case 'latinhypercube'
        % Generate a latin hypercube sample for every iteration
        srchSpc = zeros(searchSize,iterations*2);
        for i = 1:iterations
            srchSpc(:,i*2-[1 0]) = lhsdesign(searchSize,2); 
        end
end

% iterate for every signal
for i = 1:numSignals
if nargout>0 || (i == plotId) || groupAlign
    
    % Main loop: searches for optimum values for the Scale and Shift
    % factors by exhaustive search over a multiresolution grid, getting
    % finer at every iteration.

    %set to back to the user input arguments (or default)
    shft = shiftRange;
    scl  = scaleRange;

    for t = 1:iterations % increase for better resolution
        % scale and shift search space
        A = scl(1) + srchSpc(:,t*2-1) * diff(scl);
        B = shft(1) + srchSpc(:,t*2) * diff(shft);
        temp = A*corr_sig_x'+repmat(B,1,corr_sig_l);
        temp(:) = interp1(X,Y(:,i),temp(:),'pchip',0);
        [dump,imax] = max(temp*corr_sig_y);
        % save optimum
        sclOpt(i)  = A(imax);
        shftOpt(i) = B(imax);
        % readjust grid for next iteration
        scl  = sclOpt(i)  + [-0.5 0.5]*diff(scl) *reduceRangeFactor;
        shft = shftOpt(i) + [-0.5 0.5]*diff(shft)*reduceRangeFactor;
    end

end % if nargout>0 || (i == plotId)    
end % for i = 1:numSignals 

if groupAlign && numSignals>1
    shftAll =  mean(mean(sclOpt*P'+repmat(shftOpt,1,numPeaks)-...
                                   repmat(P',numSignals,1)));
    shftOpt = shftOpt - shftAll;
    P = P + shftAll;
end

% iterate for every signal
for i = 1:numSignals
if nargout>0 || (i == plotId)
    
    if (i == plotId)
        yo = Y(:,i);
    end
    
    % interpolation back to the original domain
    Y(:,i) = interp1((X-shftOpt(i))/sclOpt(i),Y(:,i),X,'pchip',NaN);
    
    if (i == plotId)
        figure
        plot(X,[yo Y(:,i)]);
        hold on
        han = plot([P P]',[min(Y(:,i));max(Y(:,i))]*ones(1,numPeaks),'r--','Tag','Xmarker');
        for j = 2:numel(han)
            setappdata(han(j),'legend_hgbehavior',0)
        end
        title(sprintf('Signal ID: %d',i));
        xlabel('Separation Units')
        ylabel('Relative Intensity')
        legend('Original signal','Aligned signal',...
               'Reference peaks')
        axis([min(X) max(X) min(Y(:,i)) max(Y(:,i))])
        grid on
        hold off
        setAllowAxesRotate(rotate3d(gcf),gca,false)
    end
    
end % if nargout>0 || (i == plotId)    
end % for i = 1:numSignals 

if nargout == 0
    clear Y
end

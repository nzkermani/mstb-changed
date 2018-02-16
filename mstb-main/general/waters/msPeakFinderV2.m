function Peaks = msPeakFinderV2( mzVals, mzInts, varargin )
%% MSPEAKFINDERV2 finds peaks in raw mass spectra
%
%   mzVals: vector containing the m/z values (num variables x 1)
%   mzInts: matrix containing the spectrum intensities (num variables x num
%           samples)
%   options:
%           - 'baseline', true/false (default true): baseline correction
%           - 'smoothing', true/false (default false): smoothing preprocess
%           - 'integrate', true/false (default false): integrate peaks
%           - 'display', true/false (default false): plot 1 random spectrum
%   if 'smoothing' is true:
%   - 'smoothMethod' (default 'sgolay'):
%      1. 'sgolay': use Savitzky-Golay filter
%          options for 'sgolay':
%          a. 'order' (default 3):        n (polinomial order for the filter)
%          b. 'frameLength' (default 11): f (frame size)
%
%      2. 'wavelet': use wavelet denoising
%          options for 'wavelet':
%          a. 'thresRule' (default 'heursure'): 'rigrsure', 'heursure', 'sqtwolog',
%                                               'minimaxi' (threshold selection rule)
%          b. 'softHard' (default 's'):         's', 'h' (soft or hard thresholding)
%          c. 'rescaling' (default 'one'):      'one' (no rescaling),
%                                               'sln' (rescaling using a single estimation of
%                                               level noise based on first-level coefficients)
%                                               'mln' (rescaling done using level-dependent
%                                               estimation of level noise)
%          d. 'level' (default 3):              n (wavelet decomposition level, must be a
%                                               positive scalar)
%          e. 'wfilter' (default 'sym8'):      'haar' (Daubechies),
%                                              'dmey' (discrete Meyer),
%                                              'sym2'...'sym8' (symlets)
%
% Examples.
% Peaks = msPeakFinderV2(mz, Y);      finds the peaks in Y given
% the m/z vector. No smoothing preprocessing is applied. Peaks are baseline
% corrected. Peaks integrals are not computed.
%
% Peaks = msPeakFinderV2(mz, Y, 'smoothing', true, 'integrate', true);
% finds the peaks in Y given the m/z vector. Spectra are smoothed using the
% Savitzky-Golay filter, with an order of 3 and a frame length of 11. Peaks
% are baseline corrected. Peak integrals are computed.
%
% Peaks = msPeakFinderV2(mz, Y, 'baseline', false, 'smoothing', true, ...
% 'smoothMethod', 'wavelet', 'thresRule', 'minimaxi', 'level', 5, ...
% 'wfilter','haar');
% finds the peaks in Y given the m/z vector. Spectra are smoothed using the
% wavelet filter, with the 'minimax' threshold rule, the threshold is soft,
% no rescaling. The wavelet decomposition is made at level 5 using the
% Haar wavelet. Peaks are not baseline corrected, and peak integrals are
% not computed.
%

options         = getVarArgin(varargin);
%strOpt          = options2Strings(options);

nVarbls         = size(mzInts, 1);
nObs            = size(mzInts, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% display options
%fprintf('------------------\n');
%fprintf('Peak finder script\n');
%fprintf('------------------\n');
%fprintf('Smoothing: %s\n', strOpt.smoothing);
%if (options.smoothing)
%    fprintf('- Smoothing method: %s\n', options.smoothMethod);
%    switch options.smoothMethod
%        case 'sgolay'
%            fprintf('- Smooth order: %d\n', options.sgolay.order);
%            fprintf('- Frame size: %d\n', options.sgolay.frameLength);
%        case 'wavelet'
%            fprintf('- Threshold rule: %s\n', options.wavelet.thresholdRule);
%            fprintf('- Soft/hard threshold: %s\n', options.wavelet.softHard);
%            fprintf('- Rescaling: %s\n', options.wavelet.rescaling);
%            fprintf('- Decomposition level: %d\n', options.wavelet.level);
%            fprintf('- Wavelet filter: %s\n', options.wavelet.wfilter);
%    end
%end
%fprintf('Baseline correction: %s\n', strOpt.baseline);
%fprintf('Peak integration: %s\n', strOpt.integrate);
%if (options.integrate)
%    fprintf('- Integral method: %s\n', options.integralMethod);
%end
%fprintf('Plots: %s\n', strOpt.display);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% do smoothing
if (options.smoothing)
    %%% perform smoothing
    switch options.smoothMethod
        case 'sgolay'
            mzInts = sgolayfilt(mzInts, options.sgolay.order, options.sgolay.frameLength);
            zeroIntIndcs = (mzInts ==0);
            mzInts(mzInts < 0) = 0; %%% sometimes the filter can give negative intensities
        case 'wavelet'
            mzInts = wden(mzInts, options.wavelet.thresholdRule, options.wavelet.softHard, ...
                options.wavelet.rescaling, options.wavelet.level, options.wavelet.wfilter);
            mzInts = reshape(mzInts, nVarbls, nObs);
            mzInts(mzInts < 0) = 0; %%% sometimes the filter can give negative intensities
    end
end

%% peak finder

nVrbls         = length(mzVals);
% nSamples       = size(mzInts, 2);

%%% fix the cuts
mzInts([1,2,nVrbls-1,nVrbls], :) = 0;

d = diff(mzInts, 1, 1); %%% compute first derivative of spectra
xderiv = [zeros(1, size(mzInts, 2)); d]; %%% put a 0 for the first element of the derivative

% adding minor perturbations & fluctuations to avoid zero derivate signals
eps = 10.^-4;
zeroderiv  = xderiv==0& zeroIntIndcs~=0;
if sum(zeroderiv)~=0
    mzInts(zeroderiv) = mzInts(zeroderiv) + cumsum(rand(length(xderiv(zeroderiv)),1))*min(abs(mzInts(zeroderiv~=0)))*eps;
    % recalculating the derivative
    mzInts([1, nVrbls], :) = 0;
    d = diff(mzInts, 1, 1); %%% compute first derivative of spectra
    xderiv = [zeros(1, size(mzInts, 2)); d]; %%% put a 0 for the first element of the derivative
end

%%% conditions for peak finding
zeroToPosSlope = (xderiv==0 & xderiv([2:nVrbls,1], :)>0);
negToPosSlope  = (xderiv<0 & xderiv([2:nVrbls,1], :)>0);
negToZeroSlope = (xderiv<0 & xderiv([2:nVrbls,1], :)==0);
posToZeroToPos = (xderiv>0 & xderiv([2:nVrbls, 1], :) == 0 & xderiv([3:nVrbls, 1, 1], :) > 0);
negToZeroToNeg =  (xderiv<0 & xderiv([2:nVrbls, 1], :)==0 & xderiv([3:nVrbls, 1, 1], :)<0);

%%% start peak
peakStartPos  = ((zeroToPosSlope) | (negToPosSlope));
if nnz(posToZeroToPos) > 0
    [badIdx, badSpec] = find(posToZeroToPos);
    nBadSpectra = length(unique(badSpec));
    for iB = 1:nBadSpectra
        badPeaks = badIdx(badSpec == iB);
        nBadPeaks = length(badPeaks);
        startIdx = find(peakStartPos(:, iB));
        for jB = 1:nBadPeaks
            peakStartPos(startIdx(find(startIdx < badPeaks(jB), 1, 'last')+1)) = false;
        end
    end
end

% mzStartPos = mzMatrix(peakStartPos);
% mzStartInt = mzInts(peakStartPos);

%%% end peak
peakEndPos    = ((negToZeroSlope) | (negToPosSlope));
if nnz(negToZeroToNeg) > 0
    peakEndPos(negToZeroToNeg) = false;
end

% Sort out the end positions if unequal from the start position quantity
jS = find(peakStartPos);
jE = find(peakEndPos);
if numel(jS) ~= numel(jE)
    
    all = unique([jS; jE]);
    
    sJ = all(1:end-1);
    eJ = all(2:end);
    
    sz = size(peakStartPos);
    peakStartPos = false(sz);
    peakEndPos = false(sz);
    peakStartPos(sJ) = true;
    peakEndPos(eJ) = true;
    
    % Now find the maximal intensities between these two points
    peakMaxPos = false(sz);
    for n = 1:numel(sJ)
        
        [~,b] = max(mzInts(sJ(n):eJ(n)));
        c = sJ(n) + b - 1;
        peakMaxPos(c,1) = true;
    end
    

else
    %%% max peak
    peakMaxPos = (xderiv >= 0 & xderiv([2:nVrbls,1], :) < 0);
    if nnz(negToZeroToNeg) > 0
        [badIdx, badSpec] = find(negToZeroToNeg);
        nBadSpectra = length(unique(badSpec));
        for iB = 1:nBadSpectra
            badPeaks = badIdx(badSpec == iB);
            nBadPeaks = length(badPeaks);
            MaxIdx = find(peakMaxPos(:, iB));
            for jB = 1:nBadPeaks
                peakMaxPos(MaxIdx(find(MaxIdx > badPeaks(jB), 1, 'first'))) = false;
            end
        end
    end
end
    

%%%%%%%%%%%%%%%%%%%%%
% REALLY IMPORTANT! %
%%%%%%%%%%%%%%%%%%%%%
%%% check if the number of peaks start, end, and max positions are the same
if  (nnz(peakStartPos) ~= nnz(peakEndPos)) %%% we just check different number of start end pos
    warning('number of peak start, end, and max positions are not the same.');
    
    %%% quick find the samples to be corrected
    StartEndPos = peakStartPos + peakEndPos;
    numPositions = sum(StartEndPos, 1);
    badSamples = find(mod(numPositions, 2) ~= 0);
    
    %%% remove the start or the end positions in excess (also the
    %%% corresponding maximum if there is)
    for j = 1:length(badSamples)
        if (nnz(peakMaxPos(:, badSamples(j))) > nnz(peakEndPos(:, badSamples(j))))
            peakMaxPos(find(peakStartPos(:, badSamples(j)), 1, 'last'):end, badSamples(j)) = false;
        end
        if (nnz(peakStartPos(:, badSamples(j))) > nnz(peakEndPos(:, badSamples(j))))
            peakStartPos(find(peakStartPos(:, badSamples(j)), 1, 'last'), j) = false;
        end
        
       if(nnz(peakStartPos(:, badSamples(j))) > nnz(peakEndPos(:, badSamples(j))))
            %%% peakMaxPos should be fixed before because it uses
            %%% peakEndPos           
            
            peakStartPos(find(peakStartPos(:, badSamples(j)), 1, 'last'), j) = false;
        elseif (nnz(peakStartPos(:, badSamples(j))) < nnz(peakEndPos(:, badSamples(j))))
            if (nnz(peakMaxPos(:, badSamples(j))) < nnz(peakEndPos(:, badSamples(j))))
                peakMaxPos(find(peakStartPos(:, badSamples(j)), 1, 'last'):end, badSamples(j)) = false;
            end
            
            peakMaxPos(1:find(peakEndPos(:, badSamples(j)), 1, 'first'), badSamples(j)) = false;
            peakEndPos(find(peakEndPos(:, badSamples(j)), 1, 'first'), badSamples(j)) = false;
        end
    end
end

%% do baseline correction
mzIntsCorr = []; %%% by default it's empty
if(options.baseline)
    %%% define a matrix with all the mz values for the interpolation
    %%% columns contain mz values
    %%% because the find command works column-wise
    mzMatrix = repmat(mzVals, 1, size(mzInts, 2));
    mzShift = [0, cumsum(repmat(mzVals(end), 1, size(mzInts, 2)-1))];
    mzMatrix = bsxfun(@plus, mzMatrix, mzShift);
    
    mzStartEndPos = mzMatrix(peakStartPos | peakEndPos);
    mzStartEndInt = mzInts(peakStartPos | peakEndPos);
    
    %%% query point for all the samples
    ps = find(peakStartPos);
    pe = find(peakEndPos);
    
    xqTmp = zeros(1, max(pe)+1);
    xqTmp(ps) = 1; %#ok<FNDSB>
    xqTmp(pe+1) = xqTmp(pe+1)-1;
    xqTmp = find(cumsum(xqTmp)); 
    
    xqTmp(xqTmp>length(mzMatrix))= []; % KV ad hoc bug fix
    xq = mzMatrix(xqTmp);
    
    %%% griddedInterpolant is much faster than interp1
    %yq = interp1(mzStartEndPos, mzStartEndInt, xq); %xq(~isnan(xq)));
    F = griddedInterpolant(mzStartEndPos, mzStartEndInt);
    yq = F(xq);
    
    %%% baseline matrix
    Baseline = zeros(size(mzInts));
    Baseline(xqTmp) = yq;
    
    %%% baseline correction
    mzIntsCorr = mzInts - Baseline;
    mzIntsCorr(mzIntsCorr < 0) = 0;
end

%% do integration (rectangular)
if (options.integrate)
    xAxis          = (1:nVarbls)';
    deltaMz        = [xAxis(2:end) - xAxis(1:end-1); 0];
    switch options.integralMethod
        case 'rectangular'
            if ~isempty(mzIntsCorr)
                totalIntegral = cumsum(mzIntsCorr);
            else
                totalIntegral = cumsum(mzInts .* deltaMz(:, ones(1, size(mzInts, 2))));
            end
            peakIntegralTmp = totalIntegral(peakEndPos) - totalIntegral(peakStartPos);
            
        case 'trapezoid'
            %%% DOUBLE CHECK THIS !!!
            if ~isempty(mzIntsCorr)
                totalIntegral = cumsum([zeros(1, size(mzIntsCorr, 2)); (mzIntsCorr(2:end, :)+mzIntsCorr(1:end-1, :))]/2 .* deltaMz(:, ones(1, size(mzIntsCorr, 2))));
            else
                totalIntegral = cumsum([zeros(1, size(mzInts, 2)); (mzInts(2:end, :)+mzInts(1:end-1, :))]/2 .* deltaMz(:, ones(1, size(mzInts, 2))));
            end
            peakIntegralTmp = totalIntegral(peakEndPos) - totalIntegral(peakStartPos);
    end
    peakIntegral = zeros(size(mzInts));
    peakIntegral(peakMaxPos) = peakIntegralTmp;
end

%% plot one spectrum for visual inspection
%%% visualize just 1 random spectrum
visSpectra = 1;
if (options.display)
    if size(mzInts, 2) > visSpectra
        visSpectra = randsample(size(mzInts, 2), visSpectra);
    end
    
    plot(mzVals, mzInts(:, visSpectra), 'b'); hold on
    if (options.baseline)
        % plot(mzVals, peakIntegral(:, visSpectra), 'color', [0.1 1 0.1]);
        %legend('raw spectrum', 'corr. spectrum (integral peaks)', 'baseline');
        %else
        peakMaxIndcs = find(peakMaxPos(:, visSpectra));
        plot(mzVals(peakMaxIndcs), mzIntsCorr(peakMaxIndcs, visSpectra), 'ro');
        plot(mzVals, mzIntsCorr(:, visSpectra), 'color', [0.1 1 0.1]);
        plot(mzVals, Baseline(:, visSpectra), ':');
        legend('raw', 'max corrected', 'corrected', 'baseline');
    end
    hold off   
    xlabel('m/z');
    ylabel('intensity');
    title(['spectrum #', num2str(visSpectra)]);
end

%% build the output
Peaks = cell(size(mzInts, 2), 1);
for p = 1:size(mzInts, 2)
    Peaks{p}.mz = mzVals;
    Peaks{p}.numPeaks = nnz(peakMaxPos(:, p));
    Peaks{p}.peakStartPos = find(peakStartPos(:, p));
    Peaks{p}.peakEndPos = find(peakEndPos(:, p));
    Peaks{p}.peakMaxPos = find(peakMaxPos(:, p));
    if (options.integrate)
        Peaks{p}.peakIntegral = peakIntegral(Peaks{p}.peakMaxPos, p);
    end
    if (options.baseline)
        Peaks{p}.baseline = Baseline(:, p);
        Peaks{p}.peakCorrectIntensity = mzInts(:, p);
    end
    Peaks{p}.mzPeaks = Peaks{p}.mz(Peaks{p}.peakMaxPos);
    if (options.integrate)
        Peaks{p}.peakInt = Peaks{p}.peakIntegral;
    else
        Peaks{p}.peakInt = Peaks{p}.peakCorrectIntensity(Peaks{p}.peakMaxPos);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [options] = getVarArgin(argsin)

%%% smoothing = true, false
%%% smoothMethod = 'sgolay', 'wavelet'
%%% baseline = true, false
%%% integrate = true, false
%%% integralMethod = 'rectangular', 'trapezoid'
%%% display = true, false

options.baseline = true;
options.smoothing  = false;
options.smoothMethod = 'sgolay';
options.integrate = true;
options.integralMethod = 'trapezoid';
options.display = false;

%%% default options for sgolay smoothing
options.sgolay.order = 3;
options.sgolay.frameLength = 11;

%%% default options for wavelet smoothing
options.wavelet.thresholdRule = 'heursure';
options.wavelet.softHard = 's';
options.wavelet.rescaling = 'one';
options.wavelet.level = 3;
options.wavelet.wfilter = 'sym8';

nArgs   = length(argsin);
for i=1:2:nArgs
    
    switch argsin{i}
        
        case 'smoothing'
            if islogical(argsin{i+1})
                options.smoothing  = argsin{i+1};
            else
                error('smooth must be true or false');
            end
            continue;
            
        case 'smoothMethod'
            switch argsin{i+1}
                case 'sgolay'
                    options.smoothMethod = 'sgolay';
                    options.smoothing = true;
                    continue;
                case 'wavelet'
                    options.smoothMethod = 'wavelet';
                    options.smoothing = true;
                    continue;
                otherwise
                    error('smoothMethod must be ''sgolay'' or ''wavelet''');
            end
            continue; %#ok<UNRCH>
            
        %%% options for sgolay smoothing    
        case 'order'
            if isscalar(argsin{i+1})
                if ~strcmp(options.smoothMethod, 'sgolay')
                    error('order option must be used with sgolay smoothing');
                else
                    options.sgolay.order = argsin{i+1};
                    options.smoothing = true;
                end
            else
                error('smoothing order must be a scalar');
            end
            continue;
            
        case 'frameLength'
            if isnumeric(argsin{i+1})
                options.frameLength = argsin{i+1};
            else
                error('frameLength order must be a scalar');
            end
            continue;
            
        %%% options for wavelet smoothing
        case 'level'
            if isnumeric(argsin{i+1}) && isscalar(argsin{i+1})
                if ~strcmp(options.smoothMethod, 'wavelet')
                    error('level option must be used with wavelet smoothing');
                else
                    options.wavelet.level = argsin{i+1};
                    options.smoothing = true;
                end
            else
                error('level must be a scalar');
            end
            
        case 'thresRule'
            if strcmp(argsin{i+1}, 'rigrsure') || ...
                    strcmp(argsin{i+1}, 'heursure') || ...
                    strcmp(argsin{i+1}, 'sqtwolog') || ...
                    strcmp(argsin{i+1}, 'minimaxi')
                if ~strcmp(options.smoothMethod, 'wavelet')
                    error('thresRule option must be used with wavelet smoothing');
                else
                    options.wavelet.thresholdRule = argsin{i+1};
                    options.smoothing = true;
                end
            else
                error('thresRule can be ''rigrsure'', ''heursure'', ''sqtwolog'', ''minimaxi''');
            end
            
        case 'rescaling'
            if strcmp(argsin{i+1}, 'one') || ...
                    strcmp(argsin{i+1}, 'sln') || ...
                    strcmp(argsin{i+1}, 'mln')
                if ~strcmp(options.smoothMethod, 'wavelet')
                    error('rescaling option must be used with wavelet smoothing');
                else
                    options.wavelet.rescaling = argsin{i+1};
                    options.smoothing = true;
                end
            else
                error('rescaling option can be ''one'', ''sln'', ''mln''');
            end
            
        case 'softHard'
            if strcmp(argsin{i+1}, 's') || ...
                    strcmp(argsin{i+1}, 'h')
                if ~strcmp(options.smoothMethod, 'wavelet')
                    error('softHard option must be used with wavelet smoothing');
                else
                    options.wavelet.softHard = argsin{i+1};
                    options.smoothing = true;
                end
            else
                error('softHard can be ''s'', ''h''');
            end
            
        case 'wfilter'
            if strcmp(argsin{i+1}, 'haar') || ...
                    strcmp(argsin{i+1}, 'dmey') || ...
                    strcmp(argsin{i+1}, 'sym2') || ...
                    strcmp(argsin{i+1}, 'sym3') || ...
                    strcmp(argsin{i+1}, 'sym4') || ...
                    strcmp(argsin{i+1}, 'sym5') || ...
                    strcmp(argsin{i+1}, 'sym6') || ...
                    strcmp(argsin{i+1}, 'sym7') || ...
                    strcmp(argsin{i+1}, 'sym8')
                if ~strcmp(options.smoothMethod, 'wavelet')
                    error('wfilter option must be used with wavelet smoothing');
                else
                    options.wavelet.wfilter = argsin{i+1};
                    options.smoothing = true;
                end
            else
                error('wfilter can be ''haar'', ''dmey'', ''sym2'', ''sym3'', ''sym4'', ''sym5'', ''sym6'', ''sym7'', ''sym8''');
            end
                    
        %%% options for peak integration
        case 'integrate'
            if islogical(argsin{i+1})
                options.integrate  = argsin{i+1};
            else
                error('integrate must be true or false');
            end
            continue;
            
        case 'integralMethod'
            switch argsin{i+1}
                case 'rectangular'
                    options.integralMethod = 'rectangular';
                    options.integrate = true;
                    continue;
                case 'trapezoid'
                    options.integralMethod = 'trapezoid';
                    options.integrate = true;
                    continue;
                otherwise
                    error('integralMethod must be ''rectangular'' or ''trapezoid''');
            end
            continue; %#ok<UNRCH>
            
        case 'display'
            if islogical(argsin{i+1})
                options.display  = argsin{i+1};
            else
                error('display must be true or false');
            end
            continue;
            
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OptString] = options2Strings(options)

if (options.smoothing)
    OptString.smoothing = 'true';
else
    OptString.smoothing = 'false';
end

if (options.display)
    OptString.display = 'true';
else
    OptString.display = 'false';
end

if (options.integrate)
    OptString.integrate = 'true';
else
    OptString.integrate = 'false';
end

if (options.baseline)
    OptString.baseline = 'true';
else
    OptString.baseline = 'false';
end

end

function [cmz,SpPeaks,mzpeaksMatched,matchedindcs] = mspmatch(mzpeaks,varargin)
% mspmatch matches mass spectrometry peaks from different spectral profiles.
%
%  [cmz pout] = mspmatch(mzpeaks) matches peak lists from multiple spectral
%  profiles(peak picked data). The algorithm identifies a common m/z vector across
%  all spectra by using a kernel density approach and then performs
%  a mathematically optimal matching of m/z values of a spectrum
%  to this common m/z vector by means of dynamic programming or nearest
%  neighbour approach

%  Input:
%         plist [1 x number of samples] is a cell vector of MS peak picked data, where each cell
%         contains peak mass-charge values and the second column has peak ion intensities
%         (peak maximum or integral values)

%  mspmatch(...,'display',value) displays the results of common mz vector
%                                estimatetion when value sets to 1 (default)
%  mspmatch(...,'ESTIMATIONMETHOD',EM) sets the method used to find the
%  common mass-charge vector (CMZ). Options are:
%
%    'histogram'    - Peak locations are clustered using a kernel density
%     (default)       estimation approach. The peak ion intensity is used
%                     as a weighting factor. The center of all the found
%                     clusters conform the CMZ vector.
%
%    'regression'   - Takes a sample of the distances between observed
%                     significant peaks and regresses the inter-peak
%                     distance to create a CMZ vector with a similar
%                     inter-element distances.
%
%
%  mspmatch(...,'matchmethod',CM) sets the method to align each peak picked
%  spectrum profile to the CMZ.
%
%    'shortest-path' - For each common peak in the CMZ vector, its
%                      counterpart in each peak list is selected using the
%                      shortest path algorithm.
%   'nearest-neighbor' - For each common peak in the CMZ vector, its
%        (default)        counterpart in each peak list is the peak that is
%                         closest to the common peak's m/z value.

%   Copyright 2006-2008 The MathWorks, Inc.
%
%
% References:
% [1] Jeffries, N. "Algorithms for alignment of mass spectrometry proteomic
%     data", Bioinformatics, 21(14):3066-3073, 2005.
%
% [2] Purvine, S. et.al. "Spectral Quality Assessment for High-Throughput
%     Tandem Mass Spectrometry Proteomics", OMICS A Journal of Integrative
%     Biology 8(3), 2004.



% Get input arguments
[opts] = getVarArgin(varargin);

% Decide if we are going to reload various things in the file?
if length(opts.progressFile) > 2 && exist(opts.progressFile,'file')
    % Reload all parts from the file now and ensure that they aren't
    % created/overwritten later on
    resume = load(opts.progressFile);
   
    % Overwrite opts...
    oldOpts = opts;
    opts = resume.opts;
end
    
% Estimate vector of common mass to charge values
if isempty(opts.cmz)
    
    % RUn the cmz determination function
    [cmz(:,1),cmz(:,2)] = getcmz(mzpeaks,...
        opts.intThr,...
        opts.mzRes,...
        opts.estimationMethod,...
        opts.display,...
        opts.cmzThreshold,...
        opts.thresholdMethod);
    
    % Scale intensities between 0 and 1
    cmz(:,2) = cmz(:,2) / max(cmz(:,2));
    drawnow;
    
    % For the msipmatch function, we determine the cmz vector and then
    % quit. Then we come back to this function to run the rest of it, i.e.
    % perform sample aligmnent on each pixel
    if nargout == 1
        return
    else
        opts.cmz = cmz;
    end
    
elseif opts.progress
    % Do nothing
else
    if sum(size(opts.cmz) == 1) == 1
        opts.cmz = [opts.cmz' ones(numel(opts.cmz),1)];
    end
end

% Memory pre-allocation
if opts.progress && exist('resume','var')
    mzpeaksMatched = resume.mzpeaksMatched;
    matchedindcs   = resume.matchedindcs;
    indStart       = resume.i + 1;
    disp(['%%%% Resuming from spectrum number = ' int2str(indStart) ' %%%%']);
else
    mzpeaksMatched = cell(size(mzpeaks));
    matchedindcs   = cell(size(mzpeaks));
    indStart       = 1;%2865;     
end
nSmpls         = length(mzpeaks);

% There are two ways to run the sample aligmnent function...
switch opts.correctionMethod
    
    % Shortest path is the dynamic programming routine that is being
    % implemented here, so make the code and comments nicely transparent
    case 'shortest-path'
      
        % rpks - this is the cmz vector and the frequencies included from
        % the determination of the cmz vector. The frequencies serve as a
        % proxy for visualising intensity, but are NOT used in the
        % determination of the optimal matching solution.
        rpks = opts.cmz;
        rpks(:,2) = 0;
        
        % Use the provided maxPeakshift (or make one...)
        if isempty(opts.maxPeakshift)
            opts.maxPeakshift = median(diff(opts.cmz(:,1)));
        end
        
        % Before we continue we need to make the function handles for the
        % determination of band and gap which depend on the values in the
        % uiTree
        opts.handBand = eval(['@(z) (' num2str(opts.maxPeakshift) '* z / 1e6)']);            
        opts.handGap  = eval(['@(z) (' num2str(opts.gapPenalty)   '* z(:,1) / 1e6)']);
        
        %opts.handBand = 0.1;
        
        % Loop through each of the samples
        for i = indStart:nSmpls
            
            % pks is the stuff to be matched to the cmz
            pks = sortrows(mzpeaks{i});
            
            % Some things are empty, so just skip these...
            if numel(pks) == 2
                % Then cannot do the samplealign function. Instead use the
                % max peak shift to determine to which cmz value the peak
                % is closest. If none, then ignore, otherwise include...
                
                % Determine to which mz values we can pair the pks(1)
                mzrange = [pks(1)-opts.maxPeakshift pks(1)+opts.maxPeakshift];
                
                chk = opts.cmz(:,1) > mzrange(1) & cmz(:,1) < mzrange(2);
                
                if sum(chk) == 0
                    % do nothing
                elseif sum(chk) == 1
                    % single good match
                    j = find(chk == 1);
                    mzpeaksMatched{i} = [rpks(j,1) pks(2)];
                    matchedindcs{i}   = [j,1];
                else
                    [~,dind] = min(abs(cmz(chk) - pks(1)));
                    j = find(chk == 1);
                    j = j(dind);
                    mzpeaksMatched{i} = [rpks(j,1) pks(2)];
                    matchedindcs{i}   = [j,1];
                end
                                
            elseif ~isempty(pks)
                % b is an n row by 2 column matrix, with mz values in the
                % first column and intensities in the second. The
                % intensities are NOT used in the calculation,
                % following the modification of the original peak matching
                % function (although they ARE passed to the function).
                % The distance caluclation is now simply the absolute
                % distance between mz values.
                
                % Make it a 'full' matrix
                if issparse(pks)
                    pks = full(pks);
                end
                
                % Scale the intensity to a maximum of 1, which is strictly
                % unnecessary as intensities are not used in the distance
                % calculation
                b = double(bsxfun(@times,pks,1./[1 max(pks(:,2))]));
                
%                 if size(b,1) > 100
%                     opts.boolSA = true;
%                 else
%                     opts.boolSA = false;
%                 end
                
                % Sometimes the function fails, perhaps because there are
                % no peaks in the function and stuff, but if we catch it we
                % can it to prevent the whole thing going shitty
                try                    
                    [j,k] = samplealign2(...
                        double(rpks),...
                        b,...
                        'Band',opts.handBand,...
                        'Gap',opts.handGap,...
                        'Distance',opts.handDist,...
                        'Quantile',[],...
                        'SHOWCONSTRAINTS',opts.boolSC,...
                        'SHOWNETWORK',opts.boolSN,...
                        'SHOWALIGNMENT',opts.boolSA);
                    
                    % These are the things to be returned to help with aligning
                    % the 
                    mzpeaksMatched{i} = [rpks(j,1) pks(k,2)];
                    matchedindcs{i} = [j,k];
                    
%                     if opts.boolSA
%                         i
%                         i;
%                     end
                    
                catch error
                    error
                end
                               
            end
                        
            % Status report so we know what is happening...
            if mod(i,100) == 0
                disp([int2str(i) '/' int2str(nSmpls)]);
            end
            
            % Because Imperial college computer policy is so stupid, there
            % are not enough waking hours of this computer to complete the
            % alignment of large 3D files. Thus I have decided to implement
            % a policy of saving periodically such that we can restore up
            % to the last save point. This requires you to specify in the
            % varargin 'progress' as true, which is false by default.
            if mod(i,5000) == 0 && opts.progress
                % Save after every thousand^th alignment...
                % Need to save: opts,mzPeaksMatched,matchedindcs?,i
                tmpSav = [opts.progressPath datestr(now,'yymmdd-HHMM') '.mat'];
                disp(tmpSav);
                
                save(tmpSav,'opts','mzpeaksMatched','i','matchedindcs');
            end
        end
        
    case 'nearest-neighbor'
        % This isn't implemented well in this version
        [pks,mzpeaksMatched] = doNearestNeighbour(nSmpls,opts.cmz,mzpeaks,mzpeaksMatched);
end

% Run the alignment function...
SpPeaks = msgetpalign(opts.cmz,mzpeaksMatched);
cmz = opts.cmz;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts] = getVarArgin(argsin)
% Read and parse all of the input arguments, and make someothers by default
% if they ahve not been provided...

% CMZ
opts.estimationMethod   = 'histogram';
opts.display            = 0;
opts.cmz                = [];
opts.cmzThreshold       = 0; % 0.1?
opts.mzRes              = [];
opts.thresholdMethod    = 'normal';

% ALIGN
opts.correctionMethod   = 'shortest-path';
opts.maxPeakshift       = [];
opts.gapPenalty         = [];
opts.intThr             = 0;   
opts.handDist           = @(R,S) abs(sum((R-S),2));

% DISPLAY
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

% SAVING PROGRESS
opts.progress     = false;
opts.progressPath = '';
opts.progressFile = '';

nArgs = length(argsin);
for i=1:2:nArgs
    if strcmpi('mzres',argsin{i})
        opts.mzRes   = argsin{i+1};
    elseif strcmpi('display',argsin{i})
        opts.display = argsin{i+1};
    elseif strcmpi('correctionmethod',argsin{i})
        opts.correctionMethod = argsin{i+1};
    elseif strcmpi('estimationmethod',argsin{i})
        opts.estimationMethod = argsin{i+1};
    elseif strcmpi('cmz',argsin{i})
        opts.cmz = argsin{i+1};
    elseif strcmpi('thr',argsin{i})
        opts.intThr = argsin{i+1};
    elseif strcmpi('maxpeakShift',argsin{i})
        opts.maxPeakshift = argsin{i+1};
    elseif strcmpi('gappenalty',argsin{i})
        opts.gapPenalty = argsin{i+1};
    elseif strcmpi('cmzthreshold',argsin{i})
        opts.cmzThreshold = argsin{i+1};
    elseif strcmpi('thresholdmethod',argsin{i})
        opts.thresholdMethod = argsin{i+1};
    elseif strcmpi('progress',argsin{i})
        
        % Need to see if the file provided is an empty folder or a file? If
        % a folder, then it is the place to save. If a file, then load it
        % up and then continue saving in the same directory...
        if exist(argsin{i+1},'dir')
            opts.progressPath = argsin{i+1};
            opts.progress = true;
        elseif exist(argsin{i+1},'file')
            opts.progressFile = argsin{i+1};
            opts.progress = true;
        else
            % Fail!
        end
        
        if opts.progress            
            % Make the folder if it is specified...and nonexistent
            if length(opts.progressPath) > 2 && ~exist(opts.progressPath,'dir')
                mkdir(opts.progressPath);
            end
            
            % If we provide a mat file, determine the root folder...
            if length(opts.progressFile) > 2
                slsh = strfind(opts.progressFile,filesep);
                opts.progressPath = opts.progressFile(1:slsh(end));
            end            
        end
        
        
    else
        warning('INVALID PARAMETER OPTION');
    end
end

% Choose to display some of the graphs from dynamic programming
if opts.display == 1
    opts.boolSA = false; 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cmz,cmzJSM] = getcmz(mzPeaks,intThr,mzres,estimationMethod,...
    display,cmzThreshold,thresholdMethod)
% Estimate a common mz vector using kernel density or regression methods

% put all P into a single array
P = full(cell2mat(mzPeaks));

% limits of the MZ vector
mzmin   = min(P(:,1));
mzmax   = max(P(:,1));
mz2idx = @(x,r) round((x - mzmin) / (mzmax-mzmin) * (r-1) + .5);
idx2mz = @(x,r) (x-.5) / (r-1) * (mzmax-mzmin) + mzmin;

switch estimationMethod
    case 'histogram'
        % Count all peaks above threshold and make a coarse histogram just below
        % the resolution of the data
        if isempty(mzres)
            mzres = 0.001:0.005:1; index = 1; f = zeros(1,length(mzres));
            for i = mzres
                f(index) = getoptmzres(i,mzmax,mzmin,P,length(mzPeaks),intThr);
                index = index + 1;
            end
            optval = min(f); 
            mzres  = mzres(max(find(f==optval)));
        end
        res_mz   = round((mzmax-mzmin)/mzres)+1;
        mzva     = accumarray(mz2idx(P(:,1),res_mz),P(:,2)>intThr,[res_mz 1]);
        maxcount = max(mzva);

        if display
            mzva_bak = mzva;
        end
        
        % This defines the resolution: increase it to improve peak 
        % centroid accuracy
        K            = mzres/0.0001;        
        
        % Can decide on what to do in order to smooth / threshold the mzva
        % vector. So for thousands of 3D pixels being aligned itno a single
        % cmz vector, we will apply this method. Other methods of matching
        % around 100 samples will have a differnt method, perhaps the
        % existing one...
        switch thresholdMethod
            
            case 'jamesBaseline'
                % This approach determines the baseline value of the mzva
                % vector and subtracts it.
                smThr = mzvaSmooth(mzva)' * 3;
                
                mzva(mzva < smThr) = 0;
                mzva_bak2 = mzva;
                
            otherwise
                % The standard approach, uses simple thresholding of values
                % below 0 intensity, or the cmzThreshold as specifed in the
                % input/user parameters.
                if isempty(cmzThreshold)
                    smThr = 0;
                else
                    smThr = cmzThreshold;
                end
                
                % % Now we determine the smThr in the smoothed mzva variable
                % smThr = median(mzva(smInds))
                % 
                % % Now multiply this by the multiplier of cmzThreshold, which if
                % % blank will be set to 2
                % if isempty(cmzThreshold)
                %     smThr = smThr * 2;
                % else
                %     smThr = size(mzPeaks,1) * cmzThreshold;
                % end        
                %
                % % Here we ditch the values below the threshold...
                % mzva(mzva < smThr) = 0;

        end
        
        
        % Here we do the smoothing of the frequency/count data        
        mzva = masmooth(idx2mz(1:res_mz,res_mz),mzva,20,'loess');
        mzva(mzva < 0) = 0;
        
        % Plot the entire series of spectra / smoothed / smoo / etc in a
        % publication-style format...
        if display
            fig.fig = figure('Units','normalized',...
                'Position',[0.54 0.54 0.36 0.36],...
                'Color',[1 1 1]);
            hold on;
            
            % Determine the mz vector for mzva...
            allmz = idx2mz(1:numel(mzva),res_mz);
            
            % Original mzva
            plot(allmz,mzva_bak,'b','LineWidth',1);
            
            % Now plot the smoothed line...which is threshold multiplied
            plot(allmz,smThr,  'r','LineWidth',3);
            plot(allmz,smThr/3,'m','LineWidth',3);
            
            % Now the smThr zeroed line
            %plot(allmz,mzva_bak2,'m','LineWidth',1);
            
            % Now the smoothed version
            %plot(allmz,mzva,'k','LineWidth',2);
            
            
            xlabel('m/z',   'FontSize',18,'FontWeight','bold','FontAngle','italic');
            ylabel('Counts','FontSize',18,'FontWeight','bold');
            axis([min(allmz) max(allmz) min(mzva) maxcount]);
            box on;
            set(gca,'FontSize',16);            
            hold off;
        end
        
        
        % Interp over the mz range to get the cmz vector
        mzva         = interp1(idx2mz(1:res_mz,res_mz),mzva,idx2mz(1:res_mz*K,res_mz*K),'pchip');
        
        %figure; plot(mzva);
        cmz          = idx2mz((find(mzva(3:end)<mzva(2:end-1) & mzva(2:end-1)>mzva(1:end-2))+1),res_mz*K);
        cmzJSM       = mzva(find(mzva(3:end)<mzva(2:end-1) & mzva(2:end-1)>mzva(1:end-2))+1);
        
                
        % correct cmz for long empty spaces
        dcmz = diff(cmz);
        h    = max(0,round(dcmz./median(dcmz(:),1))-1);        
        cmz2  = [cell2mat(arrayfun(@(q) cmz(q)   +(0:h(q))*dcmz(q)/(h(q)+1),1:numel(h),'Uniform',false)) cmz(end)];
        
        % cmz is still filled with variables that are not present after the
        % smoothing - thus we may have 7000 variables remaining but these
        % are blank - time to remove them... as the above-defined threshold
        % has removed them, we only need to find the cmz values equal to
        % zero and delete those...
        
        
        if display
            msimage.mainfig.h = figure('Units','normalized',...
                'Position',[0.54 0.54 0.36 0.36],...
                'Color',[1 1 1]);
            
            hl3 = plot(cmz(ceil(1/3:1/3:numel(cmz))),...
                repmat([min(mzva) maxcount NaN],1,numel(cmz)),'r:',...
                'DisplayName','Common Mass/Charge',...
                'Tag','mzmarker',...
                'LineWidth',3);
            
            hold on;
            
            hl2 = plot(idx2mz(1:res_mz*K,res_mz*K),mzva,'-',...
                'Color',[.3 .3 1],...
                'DisplayName','Smoothed Counts',...
                'Tag','mzcounts',...
                'LineWidth',3);
            
            hl1 = plot(idx2mz(1:res_mz,res_mz),mzva_bak,'.',...
                'Color',[0 .5 0],...
                'MarkerSize',20,...
                'DisplayName','Peak Counts',...
                'Tag','mzcounts',...
                'MarkerFaceColor',[0 0.5 0]);
            
            xlabel('m/z',   'FontSize',18,'FontWeight','bold','FontAngle','italic');
            ylabel('Counts','FontSize',18,'FontWeight','bold');
            %legend([hl3,hl2,hl1])
            axis([idx2mz(1,res_mz*K) idx2mz(res_mz*K,res_mz*K) min(mzva) maxcount]);
            %grid on
            box on;
            set(gca,'FontSize',16);
            hold off
        end
        
    case 'regression'
        
        % Find all peaks above threshold and which are contiguous:
        h = find(P(:,2)>intThr);
        h = h(diff(h)==1);
        Q = [P(h,1) P(h+1,1)-P(h,1)]; % [Peak loc, Peak distance]
        Q(Q(:,2)<=0,:)=[]; % remove observations from diff spectra
        
        if size(Q,1)>10000
            Q = Q(randsample(size(Q,1),10000),:);
        elseif size(Q,1)<1000
            error(message('mspmatch:TooFewSamples'))
        end
        
        if display
            Q_bak = Q;
        end
        
        % regress a smooth curve
        Q = sortrows(Q);
        Q(:,2) = bioinfoprivate.masmooth(Q(:,1),Q(:,2),.5);
        
        % remove repeated observations in mz (so we can use interp1 later)
        Q(~diff(Q(:,1)),:)=[];
        
        % create the cmz vector
        cmz = zeros(1,round((mzmax-mzmin)/min(Q(:,2))));
        cmz(1) = mzmin;
        i = 1;
        while cmz(i)<=mzmax
            cmz(i+1) = cmz(i)+interp1(Q(:,1),Q(:,2),cmz(i),'pchip','extrap');
            i = i+1;
        end
        cmz(i:end)=[];
        
        if display
            figure
            hl3 = plot(cmz(ceil(1/3:1/3:numel(cmz))),...
                repmat([min(Q_bak(:,2)) max(Q_bak(:,2)) NaN],1,numel(cmz)),'r:',...
                'DisplayName','Common Mass/Charge','Tag','mzmarker');
            hold on
            hl1 = plot(Q_bak(:,1),Q_bak(:,2),'.','Color',[0 .5 0],...
                'MarkerSize',10,'DisplayName','Observed Distance','Tag','mzdistance');
            hl2 = plot(Q(:,1),Q(:,2),'-','Color',[.3 .3 1],...
                'DisplayName','Smoothed Distance','Tag','mzdistance');
            title('Estimated CMZ Vector by Regression Method');
            xlabel('Mass/Charge (M/Z)')
            ylabel('Inter Peak Distance')
            legend([hl3,hl2,hl1])
            axis([min(Q_bak(:,1)) max(Q_bak(:,1)) min(Q_bak(:,2)) max(Q_bak(:,2))])
            grid on
            hold off
            setAllowAxesRotate(rotate3d(gcf),gca,false)
            msDataCursor(gcf)
        end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = getoptmzres(x,mzmax,mzmin,P,nSamples,intThr)
mz2idx = @(x,r) round((x - mzmin) / (mzmax-mzmin) * (r-1) + .5);
res_mz   = round((mzmax-mzmin)/x)+1;
mzva     = accumarray(mz2idx(P(:,1),res_mz),P(:,2)>intThr,[res_mz 1]);
maxcount = max(mzva);
f        = abs(maxcount-nSamples);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SpPeaks = msgetpalign(cmz,mspeaks)
% get intensities of the matched peaks

nSmpls  = length(mspeaks);
nVrbls  = size(cmz,1);
SpPeaks = zeros(nSmpls,nVrbls);
for i = 1:nSmpls   
    if mod(i,100) == 0
        disp([int2str(i) '/' int2str(nSmpls)]);
    end
    if ~isempty(mspeaks{i})
        [vals,indcs] = intersect(cmz(:,1),mspeaks{i}(:,1)');
        SpPeaks(i,indcs) = mspeaks{i}(:,2);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pks,mzpeaksMatched] = doNearestNeighbour(nSmpls,cmz,mzpeaks,mzpeaksMatched)

lutSize = 1000000;
lut = zeros(lutSize,1);
lut([1,mz2idx(cmz(1:end-1)+(cmz(2:end)-cmz(1:end-1))/2,1000000)])=1;
lut = cumsum(lut);
for i = 1:nSmpls
    pks = mzpeaks{i};
    pks(:,1) = cmz(lut(mz2idx(pks(:,1),1000000)))';
    % remove peaks assigned to the same spot
    pks(find(~diff(pks(:,1))),2) = max([pks(find(~diff(pks(:,1))),2), pks(find(~diff(pks(:,1)))+1,2)],[],2); %#ok
    pks(find(~diff(pks(:,1)))+1,:) = [];          
    mzpeaksMatched{i} = pks;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
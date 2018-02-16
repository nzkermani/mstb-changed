function [ vec ] = genSpectra( cmz,sp )
%genSpectra - given the cmz vector create a series of spectra which differ
%in terms of mz variation, number of missing peaks, and various other
%parameters, which will then be subjected to the mspmatch function...

% localShifts - each of these will be applied to all peaks in cmz. Peaks
% will differ by ±gS(n)
localShifts = [0 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5];

% ppmShifts - this is similar to above, but based on a consistent ppm
% deviation across the mz range
ppmShifts = [5];%0.1 0.2 0.5 1 2 3 4 5 6 7 8 9 10 20 30];

numSpectra = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Before we do anything we need to do some deletion of variables in cmz and
% some in the actual spectra that will be shifted and aligned etc. This is
% easier to do before the mz values get changed...
[vec] = genGaps(cmz,sp,numSpectra);

% Now start to generate the new spectra - let's start with the global
% shifts...
%[vec.loc] = genGlobal(vec.spec(:,1),'local', localShifts);
%[vec.ppm] = genGlobal(vec.spec(:,1),'ppm',   ppmShifts);
[vec.var] = genGlobal(vec.spec(:,1),'varppm',ppmShifts);

% Save the shifted spectra to a file so we can use the same ones all the
% time, once settled on the actual method that is
%save('cmzVector.mat','vec','-append');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gen] = genGaps(mz,sp,numSpectra)
% Generate some gaps in both vectors... Do this multiple times randomly so
% we create in essence a dataset of a few (e.g. 20) pixels which will
% differ slighlty but will be substatitally the same

randDel = 15; % 5 percent of variables will be randomly deleted

% The matrix match contains details of the variables to be matched. 
% Originally, the matching will be linear, as the first cmz peak will 
% match the first in the spec, but as gaps are introduced, their values 
% are set to NaN. Note that each of the numSpectra spectra will have a
% different match vector (the second column) so this will need resetting
% each time.
master = [(1:numel(mz))' (1:numel(mz))'];

% Here we define / select the variables which will be deleted in the cmz
% Some overlap will occur, but for the most part each will contain a 
% couple of peaks which the other doesn't. Thus there will be some gaps 
% introduced to test the aligment procedure. For the generation of a large
% quantity of 'spectra', we will use random numbers to delete variables to 
% form spec. But for deletion from the cmz, we'll sort on rank and axe 
% every xth peak.
[~,rnk]         = sort(sp);
indC            = 2:15:numel(rnk);
indC            = rnk(indC);
master(indC,1)  = NaN;
numNan          = sum(isnan(master(:,1)));

% Generate the cmz vector
yes1 = find(~isnan(master(:,1)) == 1);
cmz = [mz(yes1)' sp(yes1)'];

% Here we loop through for each spectrum
spec = cell(numSpectra,1);
for s = 1:numSpectra
    
    % Restore match(:,2) to parity
    match = [master(:,1) (1:numel(mz))'];

    % Here we 'delete' some of the spec variables to introduce gaps.  We'll
    % generate a string of random numbers of lenght equal to the number of
    % variables, and 'delete' those with a number above x as a way to
    % remove 100-x% of the peaks (approximately).
    indS = rand(numel(mz),1);
    indS = find(indS < randDel / 100)
    match(indS,2) = NaN;
    
    % Save the NaN-full match matrix...
    origmatch = match;

    % Now create a spectral matrix for the variables in this spectral
    % example, i.e. pixel...
    numNan = sum(isnan(match(:,2)));
    %spec = zeros(1,numel(mz)-numNan);
    
    % Create a new matrix for matching the spectrum to the cmz vector...
    newmatch = zeros(size(match));

    % Loop through each of the variables and determine if shifts in
    % variable numbers are necessary based on inserted gaps
    for n = 1:size(match,1)

        % NaN check
        chk = isnan(match(n,:));

        % Now decide how to continue...
        if sum(chk) == 0
            % No NaNs, so this is a direct match and no need to change the
            % indices of following variables. Instead just put the
            % appropriate indices in the right newmatch matrix
            newmatch(n,:) = match(n,:);
            %cmz(n) = mz(n);            
            %spec(n) = mz(n);            

        elseif sum(chk) == 2
            % Then peak deleted from both, so can just ignore? Just need to
            % subtract one from the subsequent peak indices
            match(n:end,:) = match(n:end,:) - 1;

        elseif chk(2) == 1
            % Deleted from the spec, i.e. remains in cmz. Subtract one from
            % this column, but no need to put anything in the newmatch as
            % obvs there is no match!
            match(n:end,2) = match(n:end,2) - 1;
            %cmz(n) = mz(n);

        elseif chk(1) == 1
            % Deleted from the cmz, i.e. remains in spec. Same as above.
            match(n:end,1) = match(n:end,1) - 1;
            %spec(n) = mz(n);

        else
            % Something has gone wrong if we get here 
            error('not possible');
        end
        %clc;
        %[match(1:n+3,:) newmatch(1:n+3,:) (cmz(1:n+3))' (spec(1:n+3))']
        %n;

    end
   
    % Indices of peaks in the spectrum to keep
    yes2 = find(~isnan(origmatch(:,2)) == 1);

    % Save the spectrum with its gaps and intensity values...
    spec{s}{1} = [mz(yes2)' sp(yes2)'];
    
    % What about the matching between the two?
    zerod = find(sum(bsxfun(@eq,newmatch,[0 0]),2) == 2);
    newmatch(zerod,:) = [];
    spec{s}{2} = newmatch;
          
end
   


gen.cmz  = cmz;
gen.spec = spec;
% Now need to return the thingy to the calling function
%gen.cmz     = [cmz' cmzint'];
%gen.spec    = [spec' specint'];
%gen.match   = newmatch;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = genGlobal(cmz,type,shift)
% Generate globally shifted spectra...

% Define the function handle for calculating the shift around the cmz value
switch type
    case 'local'
        handAdd = @(x,sh) x + sh; 
        handSub = @(x,sh) x - sh;
        
    case 'ppm'
        handAdd = @(x,sh) x + (sh .* x ./ 1e6);
        handSub = @(x,sh) x - (sh .* x ./ 1e6);
        
    case 'varppm'
        handAdd = @(x,sh) x + (rand(numel(x),1) .* (sh .* x ./ 1e6));
        handSub = @(x,sh) x - (rand(numel(x),1) .* (sh .* x ./ 1e6));
                
    otherwise
        error('NO otherwise');
end

% Structure in which to store all the datasets
data = struct('type',[],'shft',[],'spec',[]);

for s = 1:numel(shift)
    
    % Copy the cmz cell which can then be shifted according to shift(s)
    spec = cmz;

    numSpec = size(spec,1);
    
    for n = 1:numSpec
        
        % For each of the spectra within this dataset generate a random
        % number to determine whether to add/subtract some mass
        numM = size(spec{n}{1});
        dec = rand(1,numM(1));
        
        % Do those values less than 0.5
        ind = dec <= 0.5;
        spec{n}{1}(ind,1) = handSub(spec{n}{1}(ind,1),shift(s));
        
        % Do those more than 0.5
        ind = dec > 0.5;
        spec{n}{1}(ind,1) = handAdd(spec{n}{1}(ind,1),shift(s));
    end
    
    % Need to save this version of 'spec' according to the shift
    data(s).type = type;
    data(s).shft = shift(s);
    data(s).spec = spec;
    
end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
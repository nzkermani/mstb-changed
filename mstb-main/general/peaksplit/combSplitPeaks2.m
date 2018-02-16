%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,MZ] = combSplitPeaks2(X,MZ,ppmTol)
% Combine splitted peaks - this version should be tailored towards the
% combination of peaks that are within the mz tolerance in PPM values,
% rather than just a fixed value...


diffmz = diff(MZ);
% if nargin < 3
%     msres = min(diffmz); 
% end
% msres = msres + msres./10;

% Calculate the 5ppm difference
ppm = ppmTol * MZ / 1e6;

% Convert the msres values into a ppm value, which will necessarily be a
% vector for each individual mz value

% Find peaks that are separated from their neighbour by (at least) the 
% ppm difference
[~,mzindcs] = find(diffmz > ppm(1:end-1));

% Determine the difference between these peaks, i.e. those that are
% different than more than the ppm tolerance
diffIND = diff(mzindcs);

% Now find differences between this vector, i.e. references to mzindcs
% which are potentially split
dubmzindcs = find(diffIND > 1);

% If there are no potentially split peaks, then we return having made no
% changes to X and MZ
if isempty(dubmzindcs)
    return;
end

% Size of X
[nrows,ncols,nvrbls] = size(X); 

% Reshape X
X = reshape(X,[nrows*ncols nvrbls]);

% Gather the variables that were determined to be unique and unsplit into
% a new variable.  Do the same for the MZ vector, which can also be trimmed
X_temp = X(:,mzindcs);

% What about taking an alternative MZ vector, that uses the mean of peaks
% that have been split and then re-joined? There is scope for this...
MZ_alt = MZ(mzindcs);
MZ_wgt = MZ(mzindcs);

% Non-complementary list of variables
noncomp = zeros(numel(dubmzindcs),1);

% Loop through each split peak and assign to its neighbouring friend...
% Does this always assign the peak to the best neighbour? It seems to
% always add to the neighbour of higher m/z value?
for i = dubmzindcs
        
    % Check that the intensity patterns are complementary, i.e. they are
    % like merged barcodes with no
    % barc = X(:,mzindcs(i)+1:mzindcs(i+1)) > 0;
    % barc2 = sum(barc,2);    
    % if any(barc2 > 1)
    %     figure; imagesc(barc);
    %     figure; stem(sum(barc,2))
    %     noncomp(i,1) = 1;
    % end

    % This combines peaks that have been split into the new X intensity
    % matrix
    X_temp(:,i+1) = sum(X(:,mzindcs(i)+1:mzindcs(i+1)),2);
    
    % This is a way to calculate a more representative MZ vector.
    % Originally the first variable was used to define the MZ value;
    % however, the combination of multiple later variables is unaccounted
    % for by this method.  By taking the average of the MZ values of the
    % split peaks, a better approximation of the MZ value can be
    % determined. This can be accomplished in two ways:
    %
    % 1) simple mean of MZ values
    MZ_alt(i+1) = mean(MZ(:,mzindcs(i)+1:mzindcs(i+1)));
    %
    % 2) weighted mean according to sparsity of mz values
    pop = X(:,mzindcs(i)+1:mzindcs(i+1)) > 0;
    pop = sum(pop,1);    
    MZ_wgt(i+1) = sum(MZ(:,mzindcs(i)+1:mzindcs(i+1)) .* pop) / sum(pop);
    
    
end

% Finally trim the mz vector
MZ = MZ(mzindcs);

% Which MZ vector to use?
MZ = MZ_wgt;

% Final return X to image size
X  = reshape(X_temp,nrows,ncols,length(mzindcs));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
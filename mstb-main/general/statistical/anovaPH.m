function [ sigNam,sigVal ] = anovaPH( mz,y,anno )
%statDiffs - starting point for visualising Nicole's bacterial data
% mz    - m/z values
% y     - spectral intensities
% anno  - cell of annotations (can be species,genus,etc... though not sure
%         of the technique's validity when compared multiple distributions

% For more information, see the following page
% http://www.tc3.edu/instruct/sbrown/stat/anova1.htm

% How many groups are there?
[unqN,~,~] = unique(anno);
numG = numel(unqN);

% How many variables
numV = size(y,2);

% Make an empty cell / double matrix to store significant peaks
sigNam = cell(numV,1);
sigVal = NaN(numV,4);

% Work through each variable at a time...
for n = 1:numV
    
    % Check that the variable isn't just zeros
    if sum(y(:,n)) > 0
        
        % Calculcate one-way ANOVA to see IF there is a significant
        % difference in the means of each groups (note that it doesn't tell
        % you which is different). The stats structure is used in the
        % post-hoc analysis later on if the p-value is deemed significant
        [p,~,stats] = anova1(y(:,n),char(anno),'off');
        
        if p == 0
            p = realmin;
        end
        
        
        % This is the post-hoc analysis which uses honestly significant
        % difference (Tukey's HSD) by default.  The outputs of multcompare
        % are: c (each row denotes a comparision betweeen groups A and B,
        % which are the first two columns. The next three columns show the
        % lower confidence limit, the mean, and the higher confidence
        % limit. If the confidence limit contains zero, the difference of
        % this comparison is not significant); m (group mean values and
        % standard errors); gn (group names, note that they are different
        % from unqN created above).
        if p < 0.05
            
            % Here do the post-hoc analysis; for description see above
            [c,m,~,gn] = multcompare(stats,'display','off');
            
            % Find the largest mean value from the first column of m. This
            % is because we want to find a peak signifying a unique
            % presence rather than a unique absence
            [~,fI] = max(m(:,1));
            
            % Find all entries in c(:,[1 2]) == fI, i.e. 3.v.fI or fI.v.4
            [fx,~] = find(c(:,1) == fI | c(:,2) == fI);

            % Check confidence limits for the presence of zero
            ll = c(fx,3) <= 0;
            hh = c(fx,5) >= 0;
            
            % How many different groups? If both ll and hh are equal to 1
            % then they span over zero; if only one then zero is not
            % within their range. So diffG represents the number of groups
            % which are statistically different to the group with the
            % highest mean intensity.
            diffG = sum(sum([ll hh],2) == 1);

            % Save for output... First the species name in a cell
            sigNam{n,1} = gn{fI};
            
            % Then the important information: number | m/z value | p-value
            % from ANOVA | number of other groups to which this variable
            % differs (identity of this group is stored in sigNam).
            sigVal(n,:) = [n mz(n) p diffG];   
            
        else
            sigVal(n,:) = [n mz(n) NaN NaN];
        end
        
    end
    
    % Progress report as it takes so long!
    if mod(n,100) == 0
        disp([int2str(n) '/' int2str(numV)]);
    end 
end

end


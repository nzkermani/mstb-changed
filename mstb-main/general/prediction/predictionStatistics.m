function [tabs] = predictionStatistics(data)
%predictionStatistics - information about the numbers of classified etc
%pixels in the rat brains

% Define the pixel thresholds that we are interested in...
threshs = [0.1:0.1:1];

numD = size(data,2);
    
% Somewhere to save...
pix = zeros(6,numD,numel(threshs));

% Variable names for the table.
varNames = cell(numD,1);

% Loop through each one...
for n = 1:numD
    
    % Save the names for the table
    varNames{n,1} = strSwap(data(n).name(1:12),' ','_');
    
    % Get the MMC probabilities
    pr = data(n).prob;
    
    % Reshape
    sz = data(n).sz;
    pr = reshape(pr,[sz(1)*sz(2) size(pr,3)]);
    
    % What is the background here?
    bg = data(n).tobg;
    
    
    % Loop through each threshold
    for r = 1:numel(threshs)
        
        % Save temporary pr with threshold applied
        tmp = pr >= threshs(r);
        
        % Sum up to decide on ambiguous/unambiguous pixels
        amb = sum(tmp(bg,:),2)  > 1;
        una = sum(tmp(bg,:),2) == 1;
        non = sum(tmp(bg,:),2) == 0;
        
        % Put the bits in the table
        pix(:,n,r) = [sum(una) sum(amb) sum(non) 0 sum(~bg) numel(bg)]';
        pix(4,n,r) = sum(pix(1:3,n,r));        
    end
end

% Create the tables for each probability threshold...
tabs = struct('tab',[]);
for n = 1:numel(threshs)
        
    % Convert the classified pixels to percentages
    pct = 100 * bsxfun(@rdivide,pix(1:4,:,n),pix(4,:,n));
    
    rN = {['Classified (' sprintf('%0.1f',threshs(n)) ')'],'Ambiguous','Unclassified','Total'};
    % Create the table
    tabs(n).tab = array2table(pct,...
        'RowNames',rN,...
        'VariableNames',varNames);
    
    disp(char(10));
    disp(tabs(n).tab);
end


end


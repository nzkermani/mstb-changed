function [ str ] = rawH5globalPick(spec)
%rawH5globalPick - stuff to plot and pick etc

tt = tic;

% How many spectra?
numS = size(spec,1);

% Lets's try to resample (slight upsample ideally) the spectra to get them
% all on the same mz vector. We don't need to use this for much, but at
% least to allow peak detection.
ss = sortrows(vertcat(spec{:,1}),1);
mzRange = [ss(1,1) ss(end,1)];

sampSize = round((size(ss,1) / numS) * 1.5);
reMZ = zeros(numS,sampSize);
resamp = zeros(numS,sampSize);

for n = 1:numS
    
    [reMZ(n,:),resamp(n,:)] = msresample(spec{n,1}(:,1),spec{n}(:,2),sampSize,...
        'Range',mzRange,...
        'ShowPlot',false);
    
end
reMZ = nanmean(reMZ,1);

% Plot all the results showing original and new signals
figure; hold on;
cols = hsv(numS);
for n = 1:numS
    
    plot(spec{n,1}(:,1),spec{n,1}(:,2),'Color',cols(n,:));
    plot(reMZ,resamp(n,:),'Color',cols(n,:));
    
end

% Save the resampled data into a structure for subsequent peak picking
str.reMZ = reMZ;
str.reSp = resamp;
str.time = toc(tt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

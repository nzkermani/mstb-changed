function [ output_args ] = h5SpecAverage(data)
%h5SpecAverage - plot the average spectra for each file in the structure...

numF = size(data,2);

figure; hold on;

cols = jet(numF);

% Loop through each file...
for n = 1:numF
    
    % Histo groups?
    numG = numel(data(n).histIDNames);
    
    for r = 1:numG
        
        % Skip background
        if strcmp(data(n).histIDNames{r},'background')
            continue;
        end
        disp(data(n).histIDNames{r});
        
        % Calculate the average
        fx = strcmp(data(n).histIDs,data(n).histIDNames{r});
        
        % Take the average and TIC normalise it
        avg = nanmean(data(n).X(fx,:));
        avg = avg / sum(avg);       
        
        % Plot
        stem(data(n).mz,avg,...
            'MarkerFaceColor',cols(n,:));
        
    end
    
end
        
    


end


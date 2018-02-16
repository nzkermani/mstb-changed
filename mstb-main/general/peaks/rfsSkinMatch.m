function [ avgs ] = rfsSkinMatch(avgs,data)
%rfsSkinMatch - from the NR15 sample, match peaks...

% This is the master spectrum
mast = avgs{3,1};
numP = size(avgs{3,1},1);
[~,ord] = sort(mast(:,2),'descend');

% Match peaks in first two samples against the last one...
for n = 1:2
    
    % Extract spectrum
    sp = avgs{n,1};
    
    % Create vector of correspondance
    vec = NaN(size(sp,1),1);
    
    % And another to signify which peaks in mast have already been taken
    take = false(numP,1);
    
    % For each of the 450 peaks in avgs{3}, we match them to peaks in sp
    for r = 1:numP
        
        % Let's start with the most intense and work that way...
        i = ord(r);
        
        % Find closest to mz(i)
        [dst,fx] = min(abs(sp(:,1) - mast(i,1)));
        if numel(fx) == 1
        else
            continue;
        end
        
        % Check that distance isn't greater than [ARBITRARY THRESHOLD]
        if dst > 0.07
            %disp('Very far away');
            
            crr = corr(data(3).sp(:,i),data(n).sp(:,fx));
            %disp(num2str(crr));
            
            if crr > 0.75
                flag = true;
            else
                flag = false;
            end
            
        else
            flag = true;
        end
            
        % Only match peaks that pass the test...
        if flag
            
            % Save to vec
            vec(fx) = i;
            
            % Set sp(fx,1) to NaN to prevent it being matched again
            sp(fx,1) = NaN;
            
            % And in take too
            take(i,1) = true;
        end
        
    end
    
    % To make sure that I haven't got it wrong, we need to plot the
    % distance between matched peaks...
    mx = [avgs{n,1}(~isnan(vec),1) mast(take,1)];
    
    sum(take)
    
    avgs{n,3} = vec;
    avgs{3,3}(:,n) = take;
    
end




end


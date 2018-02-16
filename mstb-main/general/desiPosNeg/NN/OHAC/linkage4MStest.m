function [Z, Dinternalt] = linkage4MStest(a,ppm)
%% linkage4MS calculate the linkage matrix for 1D MSI data
%% imzML data object with a given resolution (default 0.001).
%    Input:  a   -  the MS image m/z ratios of a specimen 
%            ppm -  mass tolerance of measurements, ...
%                      used to merge peaks than are closer than the mass tolerance of measurements
%    Output: Z   -  the linkage matrix [rows x 3]
%% Author: Nazanin Z. Kermani, Imperial College London 2015. 
    N = size(a,1);
    % D holds distnaces between consecitive elements of a 
    D = diff(a);
    % Eppm calculate the mass diviation for ions based of intrument
    % accuracy
    Eppm = @(mz, ppm) (2*ppm*mz)/1000000; 
    % distance between ions is scaled to their corresponding ppm error
    for j=1:N-1
        D(j)= D(j)/(Eppm(a(j),ppm/2)+Eppm(a(j+1),ppm/2));
    end
    % Dummy variable
    Dtemp =D;
    % Z is the linkage matrix 
    Z = zeros(N-1,3);
    % T keeps track of clusters' indeces
    T = 1:N;
    % Dinternal holds T from the previous iteration 
    Dinternal = T;
    Dinternalt = T;
    di = -100;
    
    %% compute the linkage matrix 
    for j=1:N-1
        if(di<=1.0000001)
        % find two data points/clusters that are closest to each other,
        % update Z (linkage matrix)
        Dinternalt = Dinternal;
        [di, i] = min(D);
        Z(j,:) = [T(i); T(i+1); di];
        % merge two data points
        if(T(i)<N && T(i+1)<N)
            Dinternal([T(i) T(i+1)]) = j+N;
            if(i>1) % check if i is not the first element of D
                D(i-1) = D(i)+D(i-1);
            end
            if(i<length(D)) % check if i is not the last element of D
                D(i+1) = D(i)+D(i+1);
            end
        else
             % merge two clusters or one data point and one cluster 
            [~ ,indxT1] = find(Dinternal == T(i));
            [~ ,indxT2] = find(Dinternal == T(i+1));
            if(T(i) < N )
                Dinternal(T(i)) = j+N;
            else
               Dinternal(Dinternal==T(i))= j+N;
            end
             if(T(i+1) < N )
                Dinternal(T(i+1)) = j+N;
            else
               Dinternal(Dinternal==T(i+1))= j+N;
            end
            if(max(indxT1) < min(indxT2))
                temp1 = indxT1;
                temp2 = indxT2;
            else
                temp1 = indxT2;
                temp2 = indxT1;
            end
            
            if(i>1) % check if i is not the first element of D
                D(i-1) = D(i-1)+sum(Dtemp(temp2-1));
            end
            if(i<length(D)) % check if i is not the last element of D
                D(i+1) = D(i+1)+sum(Dtemp(temp1));
            end
        end
        % i is the index of merged data point/cluster, exclude i and
        % continue
        T(i+1) = j+N;
        D(i)=[];
        T(i)=[];
      
    end
    end
end
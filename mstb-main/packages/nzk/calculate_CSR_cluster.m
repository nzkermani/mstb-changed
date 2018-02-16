function [P] = calculate_CSR_cluster(spImage,j,lowestPresence,numSimulations,dist)
                    %find ions coordinates on an image
                    [x y ] = find(spImage(:,:,j)>0);
                    if (size(x,1)>lowestPresence)
                        coordinates = [x y];
                        %calculate equation 3.2.16 (DD)
                        M = distance_mat(coordinates);
                        mx = max(max(M));
                        M = M + mx*eye(length(M)) ; %add to diagonal to avoid zero case
                        DD = min(M);
                        m_dist = mean(DD);
                        % options for the input of clustering sugnificance program
                        opts.m_dist = m_dist; opts.bins = 20;
                        P = sum( (  dist < opts.m_dist) + .5*(  dist == opts.m_dist) )/(numSimulations + 1);
                    else
                        P = -10000 ; 
                    end
                    clearvars -except P
             
end


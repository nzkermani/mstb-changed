function [index_notNoise,  P_value] = clusterCSR(sp,bg_pixels,noise_,numPix,lowestPresence, numSimulations,levelOfSignificant)
%% optimised CSR test for MSI
% cluster the ions based on their presence 
% for each cluset make the randomness representative once
            [nRows,nCols,nVrbls] = size(sp);
            dataCube = sp;
            sp               = reshape(sp,nCols*nRows,nVrbls);
            %% extract background pixlels
            for j=1:size(sp,2)
                topCube(:,:,j)= bg_pixels.*dataCube(:,:,j);
            end
            top_sp = reshape(topCube,nCols*nRows,nVrbls);
            %% extract a boundry around image, B{1} is the fist bound
            % make a poly in the form of 
            % first row : 1 #edges
            % second row: a b
            % ................
            % end row:    a b
            % find a more accurate way for bounadry finding
            [B,L] = bwboundaries(bg_pixels,'noholes');
            My_poly = [1 size(B{1},1) ; B{1}];
            area_ = polyarea(My_poly(:,1),My_poly(:,2));
            P_value = ones(1,length(noise_));

numClusters = unique(numPix);
distCluster = zeros( numSimulations,length(numClusters));
opts.m_dist = 1; %not used
opts.bins = 20;

parfor i = 1:length(numClusters)
    distCluster(:,i) = clust_sim(My_poly,area_,numClusters(i),numSimulations,opts,0);
end

parfor i=1:length(noise_)
    dist_ion(i) = find(numClusters ==  numPix(i));
end

parfor i=1:length(noise_)
    P_value(i) = calculate_CSR_cluster(topCube,noise_(i),lowestPresence,numSimulations,distCluster(:,dist_ion(i)));
end
  index_notNoise = (find(P_value >= 0 & P_value <levelOfSignificant)); 
  end

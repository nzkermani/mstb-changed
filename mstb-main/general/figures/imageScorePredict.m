function [ test ] = imageScorePredict(test,mmc)
%imageScorePredict - predict the scores and other things for an image given
%the MMC model of it...
%
% Careful scaling required


numF = size(test,2);


for n = 1:numF
    
    % Spectrum
    try
        dx = test(n).spAl;
        sz = test(n).szAl;
    catch
        dx = test(n).sp;
        sz = test(n).sz;
    end
        
    % Log if was done
    if mmc.opts.doLog
        dx = log(dx + test(n).logOS);                
    end
    
    % Calculate the reference spectrum for this file, including only the
    % tissue object pixels
    refIdx = test(n).tobg ~= 0;
    refSp = nanmean(dx(refIdx,:),1);
    
    
    % Determine the scores
    [T,Pr,~] = calcScoresSpectra(dx,refSp,mmc.W,mmc.B);
    
    % Reshape
    newSZ = [sz(1) sz(2) size(T,2)];
    T = reshape(T,newSZ);
    Pr = reshape(Pr,newSZ);
    
    % Put in the structure
    
    % Make a figure
    figure; imagesc(imScale(T(:,:,1:3)));
    figure; imagesc(imScale(Pr(:,:,1:3)));

    test(n).score = T;
    test(n).prob  = Pr;
    
end


end


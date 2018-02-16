function [cmz_hist,mzaligned, Spaligned,  machedindeces] = dp(mzTest, SpTest)
%% dynamic programming alignment
opts.estimationMethod = 'histogram';
opts.mzRes = 0.001;
opts.display = 0;
opts.cmz=[];
opts.correctionMethod = 'shortest-path';
opts.maxPeakshift = 2.5;
opts.gapPenalty=5;
opts.intThr = 0;

counter=1;
for i=1:size(mzTest,1)
    for j=1:size(mzTest,2) 
      mzPeaks{counter} = [squeeze(mzTest(i,j,:))   squeeze(SpTest(i,j,:))];
      counter = counter+1 ;
    end
end
   opts.cmz = mspmatch(mzPeaks',...
        'estimationMethod',opts.estimationMethod,...
        'mzres',opts.mzRes,...
        'display',opts.display);
   cmz_hist =opts.cmz;
   for i=1:size(mzTest,1)
        for j=1:size(mzTest,2)  
                mzpeaks{1} = [squeeze(mzTest(i,j,:)) squeeze(SpTest(i,j,:))];
                [cmz,SpPeaks,mzpeaksMatched,matchedindcs] = mspmatch(mzpeaks',...
                        'cmz',opts.cmz,...        
                        'correctionMethod',opts.correctionMethod,...
                        'maxPeakShift',opts.maxPeakshift,...
                        'gapPenalty',opts.gapPenalty,...
                        'display',opts.display);
%                 mzaligned(i,j).align = mzTest(i,j,matchedindcs{1}(:,2));
%                 Spaligned(i,j).align = SpTest(i,j,matchedindcs{1}(:,2));
                mzaligned(i,j).mz =mzpeaksMatched{1}(:,1);
                Spaligned(i,j).sp = mzpeaksMatched{1}(:,2);
                machedindeces(i,j).sp = matchedindcs{1}(:,1);
        end
        
   end

end
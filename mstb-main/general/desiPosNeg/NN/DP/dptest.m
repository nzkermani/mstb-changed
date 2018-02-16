opts.ppmTol = 2.5;
% Define the peak matching options
opts.estimationMethod = 'histogram';
opts.display = true;
opts.handBand = eval(['@(z) (' num2str(opts.ppmTol) '* z / 1e6)']);            
opts.handGap  = eval(['@(z) (' num2str(opts.ppmTol) '* z(:,1) / 1e6)']);
opts.handDist = @(R,S) abs(sum((R-S),2));
opts.boolSC = false; 
opts.boolSA = false; 
opts.boolSN = false;

 [dataDpNew,cmz] = mzAlignment(dataDp,opts)
% First define the cmz vector, and then align everything to it

dataDp.mz = data.mz';
dataDp.sp = data.sp';

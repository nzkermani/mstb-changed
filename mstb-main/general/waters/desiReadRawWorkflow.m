function [ pp,durn ] = desiReadRawWorkflow(sp,xy2D,rp)
%desiReadRawWorkflow - go from start to finish to get a data matrix and mz
%vector from a file...

% Time profiling...
durn = zeros(5,1);

% Convert data to cell array
tt = tic;
[data] = desiReadRaw2Image(sp,xy2D);
durn(1) = toc(tt);

% Create average spectrum
tt = tic;
[op] = imzmlBin(data,0.001);
durn(2) = toc(tt);

% Use the average spectrum and resolving power to determine a filtered
% spectrum
tt = tic;
[flt] = desiRawFilter(op.mz,op.sp,rp);
durn(3) = toc(tt);

% Determine peaks
tt = tic;
[pks] = desiRawPeakPick(op.mz,op.sp,flt,rp);
durn(4) = toc(tt);

% Extract ion images
tt = tic;
[pp.mz,pp.data] = desiRawPeakExtract(pks,data);
durn(5) = toc(tt);

%figure; bar(durn);

end


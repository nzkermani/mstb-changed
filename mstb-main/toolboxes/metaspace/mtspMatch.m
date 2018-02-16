function [dpn] = mtspMatch(dpn,mts)
%mtspMatch - match the two parts from the same file. This may involve
%rotation, etc to make them match

% Need to reshape the mts.sp structure, by matching againt the size of the
% dpn.d1.sp matrix.
szd = size(dpn.d1.sp);

% Because the metaspace image is rotated, we need to switch around the
% sizes for an easier comparison...
allsz = [mts.sz; fliplr(mts.sz)];
szdiff = abs(bsxfun(@minus,allsz,[szd(2) szd(1)]));
szdiff = sum(szdiff,2);
[~,idx] = min(szdiff);

% Let's reshape...
mts.resp = reshape(mts.sp,[allsz(idx,:) size(mts.sp,2)]);

% Rotate -90 degrees
mts.resp = rot90(mts.resp,3);

% Flip left/right
mts.resp = fliplr(mts.resp);

% What are the sizes of the images now?
szm = size(mts.resp);
if numel(szm) == 2
    szm = [szm 1];
end

% Can we trim the bottom row in mts.resp? Basically take (1,1) and fill
% outwards from there...
new = zeros([szd(1:2) szm(3)]);
new(1:szm(1),1:szm(2),:) = mts.resp(1:szm(1),1:szm(2),:);

% Calculate the difference of the two images...
%df = nansum(dpn.d1.sp,3) - nansum(new,3);
%figure; imagesc(df);

% Can we create a new dpn by replacing the various bits in dpn
dpn.d1.mz = mts.mz';
dpn.d1.sp = new;
dpn.d1.annos = mts.anno;

% So potentially we save this and can load in the toolbox...
%save(['/Users/jmckenzi/Desktop/' filename '.mat'],'dpn');

end


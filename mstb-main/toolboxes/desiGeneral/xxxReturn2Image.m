function [tmpData] = xxxReturn2Image(data,index,sz)

% What is the size of the new images?
newSz = [sz(1) sz(2) size(data,2)];

% Full size empty matrix
tmpData = NaN([newSz(1)*newSz(2) newSz(3)]);

% Restore to original location
tmpData(index,:) = data;

% Reshape
tmpData = reshape(tmpData,newSz);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

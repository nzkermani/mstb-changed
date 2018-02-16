function imgBW = doMorph(imgBW,options)
% doMorth - applies morphology operators
% 
% Opening, biggest component and closure to improve initial image 
% segmentation
%
% Input: imgBW   - image in black and white
%        options.imopen applies image opening operator if it is not empty
%        options.bigcomp retains the specified number of the biggest 
%        components if is is not empty
%        options.imclose applies image closure operator if it is not empty
% Author: Kirill Veselkov, Imperial College 2014


% Image opener operator
if  ~isempty(options.imopen)
    SE     = ones(options.imopen);
    imgBW  = myimdilate(myimerode(imgBW,SE),SE);
end

% Biggest component operator
if  ~isempty(options.bigoper)
    if options.bigoper>0
        if ~isnumeric(options.bigoper)
            nComps    = str2double(options.bigoper);
        else
            nComps    = options.bigoper;
        end
        CC              = bwconncomp(imgBW);
        numPixels       = cellfun(@numel,CC.PixelIdxList);
        nComps          = min(numel(numPixels),nComps);
        [nPix,objindcs] = sort(numPixels,'descend');
        imgBW           = zeros(size(imgBW));
        for i = 1:nComps
            imgBW(CC.PixelIdxList{objindcs(i)}) = 1;
        end
    end
end

% Image closure operator
if  ~isempty(options.imclose)
     SE     = ones(options.imclose);
     imgBW  = myimerode(myimdilate(imgBW,SE),SE);
end

end
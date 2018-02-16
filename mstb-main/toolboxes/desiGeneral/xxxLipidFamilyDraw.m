function xxxLipidFamilyDraw(~,~,fig,man)
%xxxLipidFamilyDraw - make the plots to show the world

% Get the user data
dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% Check that we have annotations
if ~isfield(dpn.d1,'anno')
    return
end

% Determine the lipid family
families = get(man.family,'String');
val = get(man.family,'Value');
opts.family = families{val};

% What about other options?
opts.suppress = get(man.suppress,'Value');
opts.log = get(man.log,'Value');

% Now write a function to extract the images from a data matrix
[~,img1] = imageExtract(dpn.d1,opts);
if strcmp(dpn.mode,'dual')
    [~,img2] = imageExtract(dpn.d2,opts);    
end

% Update the axes
dpnIonImage([],[],fig.ax.ms1,img1);
if strcmp(dpn.mode,'dual')
    dpnIonImage([],[],fig.ax.ms2,img2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx,img] = imageExtract(dx,opts)
% Extract the indices of the ions of this class...

idx = strcmp(dx.anno.fam,opts.family);
idx = sum(idx,2) == 1;

% TIC norm?
tt = nansum(dx.sp,3);
sp = bsxfun(@rdivide,dx.sp,tt) * 1000;
%sp = dx.sp;

% Quickly prepare the image of these ions
img = nansum(sp(:,:,idx),3);

% Find the 95th percentile
if opts.suppress
    pc = prctile(img(:),95);
    img(img > pc) = pc;
end

% Log transform
if opts.log
    img = log(img);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

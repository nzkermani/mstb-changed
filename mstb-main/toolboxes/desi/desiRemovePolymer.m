function desiRemovePolymer(src,event,fig)
%desiRemovePolymer - remove the polymeric signal from the data matrix.
%We'll do it in a way that let's it be added back in later on...

dpn = guidata(fig.fig);
if isempty(dpn)
    return
end

% See if we have a poly field in the dpn.d1
if ~isfield(dpn.d1,'poly')
    dpn.d1.poly.remove = false;
    dpn.d1.poly.mz = [];
    dpn.d1.poly.sp = [];
    dpn.d1.poly.idx = [];
    dpn.d1.poly.icon.noPoly = iconChange(get(src,'CData'),'red');
    dpn.d1.poly.icon.normal = get(src,'CData');
end

% Now we need to determine what the polymer peaks are
if isempty(dpn.d1.poly.mz)
    
    % Run the polymer identification function
    [pidx,~] = polymerPattern(dpn.d1.mz,dpn.d1.sp);
    
    % Remove 326.378 from H2O and OCT spectra
    fx = mzFind(dpn.d1.mz,[326.378 327.380 328.384],5);
    pidx(find(fx)) = false;
       
    % Extract the polymeric part    
    dpn.d1.poly.mz = dpn.d1.mz(~pidx);
    dpn.d1.poly.sp = dpn.d1.sp(:,:,~pidx);
    dpn.d1.poly.idx = pidx;
    
end

% Now we need to do the opposite of the stat in poly.remove. So if it says
% false then we need to remove the polymer from the main matrix.

if dpn.d1.poly.remove
    
    % Now we need to add all the polymeric peaks back in...
    newMZ = [dpn.d1.mz dpn.d1.poly.mz];
    newSP = cat(3,dpn.d1.sp,dpn.d1.poly.sp);
    
    % Now just sort the stuff
    [~,srt] = sort(newMZ);
    newMZ = newMZ(srt);
    newSP = newSP(:,:,srt);
    
    dpn.d1.mz = newMZ;
    dpn.d1.sp = newSP;
    
    % Update the icon
    dpn.d1.poly.remove = false;
    set(src,'CData',dpn.d1.poly.icon.normal);
else
    
    % Now we need to remove all the polymeric peaks
    dpn.d1.mz = dpn.d1.mz(dpn.d1.poly.idx);
    dpn.d1.sp = dpn.d1.sp(:,:,dpn.d1.poly.idx);
    
    % Update the icon
    dpn.d1.poly.remove = true;
    set(src,'CData',dpn.d1.poly.icon.noPoly);
end

% Finally update the guidata
guidata(fig.fig,dpn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ico] = iconChange(ico,colour)

mask = ico < 0.5;

i1 = ico(:,:,1);
i2 = ico(:,:,2);
i3 = ico(:,:,3);

i1(mask(:,:,1)) = 0.8;
i3(mask(:,:,3)) = 0.8;

ico(:,:,1) = i1;
ico(:,:,3) = i3;

%figure; imagesc(ico);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
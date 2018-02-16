function predesi
%predesi - function for splitting up files which contain multiple sections

% Define a default path...
if ismac
    defP = '/Volumes/JSM/DB/_Archive/3D-Liver/7TopRight, 17BottomRight, 27BottomLeft, 37TopLeft/';
    defP2 = '/Users/jmckenzi/Desktop/';
else
    defP2 = 'E:\Box Sync\PI3K1415\';
    defP = 'E:\Data\Xevo\Eicosanoid\Eiccosanoid_lock mass corrected\';
end
   
% Check default paths
if ~exist(defP,'dir') && ~exist(defP2,'dir')
    defP = [pwd filesep];
elseif exist(defP,'dir')
    
elseif exist(defP2,'dir')
    defP = defP2;
end

% Draw the window
[fig] = pdDraw;

% Update the callbacks
pdCallbacks(fig,defP);


end


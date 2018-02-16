function pdIonExtract(~,~,fig)
%pdIonExtract - extract the info from the files

% Get the file info
file = fig.fSelect.UserData;

% Now determine the ions to be extracted
ions = str2double(fig.ionList.String);
ppm  = str2double(fig.ionPPM.String);


% Which kind of file?
switch lower(file.ext)
    
    case 'imzml'        
        [extr,~,totsum] = imzmlIonExtract([file.dir file.nam],ions,ppm);
        xy2D = [];
        
    case 'raw'
        [sp,~,~,xy2D] = desiReadRaw([file.dir file.nam],true);
        [extr,~,totsum] = rawImageMake(sp,xy2D,ions,ppm);        
        
    otherwise
        error('There is no otherwise');
end

% Format the data into a single image, along with the text labels
data.image = cat(3,extr,totsum);
data.xy2D = xy2D;
data.label = [fig.ionList.String; 'TIC'];
data.ppm = ppm;
data.file = file;

% Set this to be the user data of the figure
guidata(fig.fig,data);

% Now set the ions to be displayed in the new box
fig.ionDisplay.String = data.label;
fig.ionDisplay.Value  = size(data.image,3);

% Display
pdDisplay([],[],fig);

end


function [ data,idx ] = imzmlRawExtract(file)
%imzmlRawExtract - exact all data

% Source package directory
javaclasspath(deSlash([pwd '/packages/'...
    'imzMLConverter/imzMLConverter.jar']));

% Get the file's handle thing
imzML = imzMLConverter.ImzMLHandler.parseimzML(file);

% Determine image size
nColumns = imzML.getWidth();
nRows    = imzML.getHeight();

% Create a large 3D image
data = cell(nRows,nColumns);

idx = zeros(nColumns * nRows,3);
i = 0;

% Add a waitbar
wb = waitbar(0,'Extracting');

% Run through the pixels...
for y = 1:nRows%nRows:-1:1
    for x = 1:nColumns

        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            %disp(['SKIP ' int2str(y) ',' int2str(x)]);
            continue; 
        end 

        % Get the data
        xtmp = imzML.getSpectrum(x,y).getmzArray();
        ytmp = imzML.getSpectrum(x,y).getIntensityArray();

        % Save
        data{y,x} = [xtmp ytmp];
        
        i = i + 1;
        idx(i,:) = [x y numel(xtmp)];
        

    end
    waitbar(y/nRows,wb,'Extracting');
end

% Ensure that we close the waitbar
delete(wb);

end


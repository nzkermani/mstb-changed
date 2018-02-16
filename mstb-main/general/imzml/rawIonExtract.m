function [ extr,rawmz,totsum ] = rawIonExtract(sp,xy2D,ions,ppmTol)
%rawIonExtract - extract ions and things from a raw file

% Read in the file entirely
%[sp,numP,xy,xy2D] = desiReadRaw([file.dir file.nam],true);

% Determine image size
nColumns = imzML.getWidth();
nRows    = imzML.getHeight();
nIons    = numel(ions);

% Create a large 3D image
extr = NaN(nRows,nColumns,nIons);
rawmz = NaN(size(extr));
totsum = zeros(nRows,nColumns);

% Determine the ppm tolerances of the ions...
ppm = ppmTol * ions / 1e6;

% Somewhere to store the average spectrum
%mzRes = 0.01;
%mzVector = 50:mzRes:1000;
%spec = zeros(size(mzVector));

% Add a waitbar
wb = waitbar(0,'Extracting');

% Run through the pixels...
i = 0;
for y = nRows-1:-1:1
    for x = 1:nColumns

        % Increase counter
        i = i + 1;            

        % Skip empty scans
        if isempty(imzML.getSpectrum(x,y))
            continue; 
        end 

        % Get the data
        xtmp = imzML.getSpectrum(x,y).getmzArray();
        ytmp = imzML.getSpectrum(x,y).getIntensityArray();
        
        totsum(y,x) = totsum(y,x) + nansum(ytmp);

        % Now proceed on an ion-by-ion basis
        for p = 1:nIons
            
            % Find the ions within the appropriate range
            fx = xtmp > ions(p)-ppm(p) & xtmp < ions(p)+ppm(p);
            fy = ytmp .* fx;
            
            if sum(fy) == 0
                continue;
            end
            
            % If we find something, then let's add the maximum value to the
            % matrix. This isn't very advanced, as there may be two peaks
            % within this window. Can make it more advanced in the future
            [~,locmax] = max(fy);
            
            % It would be useful to gather the peak shape and take an
            % average spectrum
            %mask = xtmp >= mzVector(1) & xtmp <= mzVector(end);
            %tmp = round((xtmp(mask)-mzVector(1)) ./ mzRes) + 1;
            %spec(tmp) = spec(tmp) + ytmp(mask)';
            
            %spec = spec + tmp;
            
            % Save to the datacube
            extr(y,x,p)  = ytmp(locmax);
            rawmz(y,x,p) = xtmp(locmax);
            
        end

    end
    waitbar(y/nRows,wb,'Extracting');
end

% Ensure that we close the waitbar
delete(wb);

% Trim spec and return the mz vector in it too
%mask = spec > 0;
%spec = [mzVector(mask)' spec(mask)'];

end


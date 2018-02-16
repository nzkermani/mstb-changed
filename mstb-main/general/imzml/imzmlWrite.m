function imzmlWrite(xml,data,opdir)
%imzmlWrite - create a new imzml and idb file

wb = waitbar(0,'Exporting');

%opdir = '/Users/jmckenzi/Desktop/imzml/';

% Determine the old binary file name (for UUID)
oldB = [xml(1:end-5) 'ibd'];

% Strip out the folder from the xml file name
sl = strfind(xml,filesep);
pth = xml(1:sl(end));
xml = xml(sl(end)+1:end);

% Create a slightly modified filename
newF = [xml(1:end-6) '-RECAL' '.imzML'];
disp(newF);

newB = [newF(1:end-5) 'ibd'];
disp(newB);

% Copy the imzml file to the new folder
copyfile([pth xml],[opdir newF]);

% Need to get the UUID, probs easier from the original binary file
fid = fopen(oldB,'r');
uuid = fread(fid,16,'uint8');

% Create new binary file
fid = fopen([opdir newB],'w');

% Now let's start to write the binary... start with uuid
fwrite(fid,uuid,'uint8');

% Now the actual data...
for p = 1:size(data,1)    
    
    for q = 1:size(data,2)
            
        fwrite(fid,data{p,q}(:,1),'double');
        fwrite(fid,data{p,q}(:,2),'single');        
        
    end 
    
    waitbar(p/size(data,1));    
    
end

fclose(fid);
delete(wb);

end


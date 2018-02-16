function [ output_args ] = dbWriteCSV(files,dir)
%dbWriteCSV - write out a CSV / txt file for metadata and stuff

fid = fopen(dir,'w');

numF = size(files,1);

fprintf(fid,'%s\n','File name');

for n = 1:numF
    
    fprintf(fid,'%s\n',files{n,1});
    
end

fclose(fid);


end


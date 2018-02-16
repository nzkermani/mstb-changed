function [binNew]= convertToBin(dims, binStruct,sorted_mz_indx,sorted_mz_value,TempSp,method)
% map sorted values back to the image
dbstop if error
method = 'restrict';
% Initialise the bin that holds the alignment results
field1 = 'mz';  value1 = zeros(1,1);
field2 = 'spectra';  value2 =zeros(1,1);
field3 = 'intensity';  value3=zeros(1,1);
field4 = 'centroidspectra';  value4 =zeros(1,1);
field5 = 'centroidmz';  value5 =zeros(1,1);

binNew = repmat(struct(field1,value1,field2,value2, field3,...
    value3,field4,value4,field5,value5),length(binStruct),1 );
tic;
h = waitbar(0,'Please wait...');
steps = length(binStruct);
counter=1;
binCounter = 1;
for i =1 :length(binStruct)
        binTemp = binStruct{1,i};
        MZ_Loc = sorted_mz_indx(counter:(counter+length(binTemp)-1));
        MZ = sorted_mz_value(counter:(counter+length(binTemp)-1));
        int = TempSp(counter:(counter+length(binTemp)-1));
        counter = counter+length(binTemp);
        j=1;
        while( j <= length(binTemp))
            indx = find(binTemp == binTemp(j));
            j = j+length(indx);
            binNew(binCounter).spectra = MZ_Loc(indx);
            binNew(binCounter).mz = MZ(indx);
            binNew(binCounter).intensity = int(indx);
%             for k = 1:length(temp1)
%                 binNew(binCounter).mz(k) = mz(temp1(k), temp2(k), temp3(k));
%                 binNew(binCounter).intensity(k) = sp(temp1(k), temp2(k), temp3(k));
%             end    
     if(strcmp(method , 'restrict'))
         
         
        [temp1, temp2, temp3] = ind2sub(dims, MZ_Loc(indx));
        if(length(temp1)>1)
         [c,ia,ic]=unique([temp1'; temp2']','rows');
        if(max(ic)<length(temp1))
            % find duplicates
            for k = 1:length(ia)
                duplicatedIndex = find(ic==ia(k));
                if(length(duplicatedIndex)>0)
                    [~, indMaxIntensity] = max(binNew(binCounter).intensity(duplicatedIndex));
                    binNew(binCounter).mz(setdiff(duplicatedIndex, indMaxIntensity))= -1;
                end
            end
        end
        end
        indExclude = find(binNew(binCounter).mz == -1);
        binNew(binCounter).mz(indExclude)=[];
        binNew(binCounter).intensity(indExclude)=[];
        binNew(binCounter).spectra(indExclude)=[];
     end
    binNew(binCounter).centroidmz = centroid( binNew(binCounter).mz, binNew(binCounter).intensity);
    binNew(binCounter).centroidspectra =  centroid( binNew(binCounter).intensity , binNew(binCounter).mz);
clear temp1 temp2 temp3
   binCounter = binCounter+1;
        end
    waitbar(i / steps);
end
toc;
close(h) 
end



function bin=binStruct2bin(binStruct)
counter = 1;
for i = 1:length(binStruct)
    for j = 1:length(binStruct{1,i})
        counter = counter+1;
    end
end
h = waitbar(0,'Prepare the aligned spectra...');
dummy = 1;
% Initialise the bin that holds the alignment results
field1 = 'mz';  value1 = zeros(1,10);
field2 = 'spectra';  value2 =zeros(1,10);
field3 = 'intensity';  value3=zeros(1,10);
field4 = 'centroidspectra';  value4 =zeros(1,1);
field5 = 'centroidmz';  value5 =zeros(1,1);
bin = repmat(struct(field1,value1,field2,value2, field3,value3,field4,value4,field5,value5),counter,1 );

for i = 1:length(binStruct)
    for j = 1:length(binStruct{1,i})
        dummy = dummy+1;
        bin(counter) = binStruct{1,i}(j);
        waitbar(dummy / counter);
    end
end
close(h)
end
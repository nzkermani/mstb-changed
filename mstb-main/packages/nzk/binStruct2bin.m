function bin=binStruct2bin(binStruct)
    counter = 1;
    for i = 1:length(binStruct)
        for j = 1:length(binStruct{1,i})
            bin(counter) = binStruct{1,i}(j);
            counter = counter+1;
        end
    end
end
%get bin and make 2d matrix
function [MZnew, SPnew] = bin2image(dims, bin, option)

nMz = numel(bin);
if(option == 'withinImage')
    SPnew = zeros(dims);
    MZnew = zeros(nMz,1);
    for i = 1:nMz
        MZnew(i,1) = bin(i).centroidmz;
        [temp1 temp2 temp3] = ind2sub(dims , bin(i).spectra);
        for (l = 1:length(bin(i).spectra))
            SPnew(temp1(l), temp2(l), i) = bin(i).intensity(l);
            
        end
    end
else
    
    SPnew = zeros(dims);
    MZnew = zeros(nMz,1);
    for i = 1:nMz
        MZnew(i,1) = bin(i).centroidmz;
        [temp1 temp2] = ind2sub(dims , bin(i).spectra);
        for (l = 1:length(bin(i).spectra))
            SPnew(temp1(l), temp2(l)) = bin(i).centroidspectra;
        end
        clear temp1 temp2
    end
end
end
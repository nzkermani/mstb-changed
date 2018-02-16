%get bin and make 2d matrix
function [MZnew, SPnew] = bin2image(mzTest, SpTest, bin, option)
  
    nMz = numel(bin);

    SPnew = zeros(size(SpTest));
    MZnew = zeros(size(SpTest));

    if(option == 'withinImage')
        SPnew = zeros(size(SpTest,1),size(SpTest,2),nMz);
        MZnew = zeros(nMz,1);
        for i = 1:nMz
             MZnew(i,1) = bin(i).centroidmz;
            [temp1 temp2 temp3] = ind2sub(size(mzTest) , bin(i).spectra);
           for (l = 1:length(bin(i).spectra))
                SPnew(temp1(l), temp2(l), i) = SpTest(temp1(l), temp2(l),temp3(l));
               
           end
        end
    else
        for i = 1:nMz
                [temp1 temp2] = ind2sub(size(mzTest) , bin(i).spectra);
               for (l = 1:length(bin(i).spectra))
                    SPnew(temp1(l), temp2(l)) = bin(i).centroidspectra;
                    MZnew(temp1(l), temp2(l)) = bin(i).centroidmz;
               end
        end
    end
end
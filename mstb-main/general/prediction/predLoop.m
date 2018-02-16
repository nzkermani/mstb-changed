function [ output_args ] = predLoop(dbdata,opts,dbP,PixelClassParam)
%predLoop - loop through each file in turn

for n = 1:size(dbdata,2)
        
    try
        [~] = dbPixelClassificationHoldOut(dbdata,opts,dbP,PixelClassParam,n);
    catch err
        err
        disp([int2str(n) 'fail']);
    end
        
    
end

end


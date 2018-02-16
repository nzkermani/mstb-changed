function closeOtherWindows(src,event)
% Close all windows other than the one that calls it...

% Is this necessary?
if rand(1) > 0.25
    try
        %[y,Fs,NBITS] = wavread('/general/desiPosNeg/icons/duck-quack4.wav');
        [y,Fs] = audioread('/general/desiPosNeg/icons/duck-quack4.wav');
        
        sound(y/4,Fs);
        
    catch
    end
end


% Handle of current figure
curF = gcf;

% Handles of all figures
allF = findobj(0,'Type','figure');

% Trim list...
delF = setdiff(allF,curF);

% And close
close(delF);

end

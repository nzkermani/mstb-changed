function [ name ] = getUser()
%getUser - simple function to return the name of the computer. Nothing
%sinister in that...
    
try
    dnam = deblank(evalc('!whoami'));    
    name = charPurge(dnam,'-');
catch
    name = 'anonymous';
end

end


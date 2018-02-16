function [logged] = logData(x)
%% logData log transforms the MS image intensity matrix 
%%    Input:  x  - Intensity matrix [rows x columns x mz]

%    Output: logged  - Intensity matrix [rows x columns x mz]
%      
%% Author: Nazanin Z. Kermani, , Imperial College London 2016.

% values smaller than 10 with be nagative if the offset is not added 
offset = min(x(x > 0));
logged = log(x + offset);
minLogged = min(logged);
if(minLogged<0)
    logged = logged - minLogged;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
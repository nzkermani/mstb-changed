function [ data ] = mass2mat( mz )
%mass2mat - convert a series of masses to columns based on thousands,
%hundreds, tens, ones, tenths, hundredths, thousandths etc...

% Round to lose extraneous decimals
mz = round(mz * 10000) / 10000;

% Define the values...
vals = [1000 100 10 1 0.1 0.01 0.001 0.0001];

% New matrix
data = zeros(numel(mz),numel(vals));

for n = 1:numel(vals)
    
    fx = mz >= vals(n);
    
    tmp = mz(fx) - mod(mz(fx),vals(n));
    
    tmp = tmp / vals(n);
    
    % Add to the matrix
    data(fx,n) = tmp;
    
    % Subtract from the mz vector
    mz(fx) = mz(fx) - (tmp * vals(n));
    
end
    
end


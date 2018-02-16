function [ op ] = imzmlBin(data,res)
%imzmlBin - generate a profile-like average spectrum from an imzML file in
%cell array format, from imzmlRawExtract

% Define an m/z vector, probably from 1 to 2000 at most
mz = 1:res:1500;
sp = zeros(size(mz));
fq = zeros(size(mz));

% Loop through
for p = 1:size(data,1)
    
    for q = 1:size(data,2)
        
        if isempty(data{p,q})
            continue;
        end
        
        tmp = data{p,q}(:,1);
        tmp = tmp - mz(1);
        
        tmp = round(tmp / res) + 1;
        
        % Check that no bigger than vector
        fx = tmp <= numel(sp) & tmp > 0;
        
        % Add in...
        %tmp2 = data{p,q}(:,2);
        sp(tmp(fx)) = sp(tmp(fx)) + data{p,q}(fx,2)';
        fq(tmp(fx)) = fq(tmp(fx)) + 1;
        
    end
    
end

op.mz = mz;
op.sp = sp;
op.fq = fq;

end


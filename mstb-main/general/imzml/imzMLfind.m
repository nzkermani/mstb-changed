function [ all ] = imzMLfind
%imzMLfind - quick function to find imzML files in a given location

extn = {'centroid.imzML','centroid.ibd',...
    'profile.imzML','profile.ibd',...
    'ndpi','jpg'};

target = '/Users/jmckenzi/Desktop/ML-ICL-Colo'

path = ['x:' genpath(target)];

sep = strfind(path,':')

for n = 2:numel(sep)-1
    
    % This is the subfolder
    tmp = path(sep(n)+1:sep(n+1)-1);
    disp(tmp);
    
    % Find all imzML & ibd & ndpi files
    dd = dir(tmp)
    numF = size(dd,1);        
    
    res = cell(1,numel(extn)+1);
    res{1,1} = tmp;

    for r = 1:numF
        
        
        % Find bit of extn in name...
        for x = 1:numel(extn)
            
            chk = strfind(dd(r).name,extn{x})
            
            if ~isempty(chk) % then we have a match
                
                res{1,x+1} = dd(r).name;
                
            end
            
                
        end
    end
       
    if ~exist('all','var')
        all = res;
    else
        all = cat(1,all,res);
    end
end
all


end


function [ head,res ] = mtspCSVimport(fname)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Just the default A30
if nargin == 0
    fpath = '/Users/jmckenzi/Dropbox/Imperial/Projects/Metaspace/ICL/Oesphageal/';
    fname = [fpath 'sm_results-6.csv'];
end

% Storage
res = cell(1000,12);

% Problem is that it isn't a full CSV with only numeric data
fid = fopen(fname);
i = 0;
while ~feof(fid)
    
    % Read the line
    tmp = [fgetl(fid) ','];
    
    if strcmp(tmp(1),'#')
        continue;
    end
    
    % Increase counter if we pass the first stage
    i = i + 1;
        
    % Find the commas
    com = [0 strfind(tmp,',')];
    
    for n = 1:12
        a = com(n) + 1;
        b = com(n+1)-1;
        
        if n >=6 && n <= 11 && i > 1
            res{i,n} = str2double(tmp(a:b));
        elseif n == 5
            res{i,n} = tmp(a+2:b-1);
        elseif n == 12
            res{i,n} = tmp(a:end-1);
        else
            res{i,n} = tmp(a:b);
        end
    end
    
    
    
        
    % Format the compound names
%     if i > 1
%         
%         % Also we need to combine the ions that have formed - because if we
%         % see +Cl and -H of the same molecule then we should treat them as
%         % separate entities rather than boshing them together
%         res(i,4) = {[res{i,3} res{i,4}]};
%         
%         % String of compound names
%         comp = tmp(com(11)+1:end);
%         
%         % Determine if there are " " tags, which suggests multiple isomeric
%         % structures
%         spm = strfind(comp,'"');
%         if numel(spm) == 2
%             
%             res{i,11} = regexp(comp(2:end-1),', ','split')';
%             
%         else
%             res{i,11} = {comp};
%         end
%         
%         
%     end
    
end

fclose(fid);

head = res(1,:);
res = res(2:i,:);

% DOn't need to determine the charge state as it isn't really important...
% But it is for the network analysis which comes later / if it gets
% resurrected.
return

% Replace ion forms of "M-H" etc with just the simple adduct part



% I want to add a column which determines the charge status of this file.
% Can't always rely on the filename, but the first ion will tell us...
posI = {'+H','+Na','+K'};
negI = {'-H','+Cl'};

pcmp1 = ~cellfun(@isempty,strfind(res(:,4),posI{1}));
pcmp2 = ~cellfun(@isempty,strfind(res(:,4),posI{2}));
pcmp3 = ~cellfun(@isempty,strfind(res(:,4),posI{3}));
pcmp = pcmp1 + pcmp2 + pcmp3;

ncmp1 = ~cellfun(@isempty,strfind(res(:,4),negI{1}));
ncmp2 = ~cellfun(@isempty,strfind(res(:,4),negI{2}));
ncmp = ncmp1 + ncmp2;

if sum(pcmp)/numel(pcmp) > 0.75
    charge = repmat({'Positive'},[size(res,1) 1]);
elseif sum(ncmp)/numel(ncmp) > 0.75
    charge = repmat({'Negative'},[size(res,1) 1]);
else
    error('Indeterminate charge state');
end
res(:,end) = charge;



end


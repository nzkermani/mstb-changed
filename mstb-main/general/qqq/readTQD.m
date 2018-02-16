function [all] = readTQD(fP,fN)
%readTQD - import a CSV from the TQ-DESI output software

if nargin == 0
    % Then we need to ask to load a file...
    [fN,fP,fI] = uigetfile('*.csv','Select CSV file','Selecta!',...
        'MultiSelect','on');
    
    if fI == 1
        path = [fP fN];
    else
        error('No file selected');
    end
end

% Run the function for multiple instances...
numF = numel(fN)
all = struct('tqd',[],'img',[],'sz',[]);
for n = 1:numF
    
    % Read in from the CSV
    [all(n).tqd] = readCSV(fP,fN{n});
    
    % Format into an image
    [all(n).img] = procTQD(all(n).tqd);
    all(n).sz = size(all(n).img);

end

% Ensure that the images are approximately the same size...
allsz = vertcat(all.sz);
minsz = min(allsz,[],1);
for n = 1:numF
    
    % Determine new image size...
    newsz = [min([allsz(n,1) minsz(1)]) min([allsz(n,2) minsz(2)]) allsz(n,3)];
    
    all(n).img = all(n).img(1:newsz(1),1:newsz(2),:);
    
end

return

% Now plot the thing!
figTQD(all);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tqd] = readCSV(fP,fN)

% Structure to store stuff
tqd = struct('name',[],'mat1',[],'mat2',[],'head',[],'data',[]);

path = [fP fN];
tqd.name = fN;

% Some counters
count0 = 0;
count1 = 0;
last = 0;


% Read in the file
fid = fopen(path,'r');
while ~feof(fid)
    
    % Read in each line
    line = fgetl(fid);
    count0 = count0 + 1;
    
    if count0 > 1
        try
            tmp = textscan(line,'%f');
        end
        
        % Now decide what to do with it
        if count0 == 2
            tqd.mat1 = tmp{1};
        
        elseif count0 == 3
            tqd.mat2 = tmp{1};
            
        elseif count0 == 4
            tqd.head = tmp{1};
            data = zeros(100000,numel(tqd.head)+3);

        else
            count1 = count1 + 1;
            data(count1,:) = tmp{1}';
        end
            
    end
    
end
   
% Trim the matrix considerably...
tqd.data = data(1:count1,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

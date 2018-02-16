function h5HistoReannotation
%h5HistoReannotation - change histIDs of tissue types according to an Excel
%spreadsheet. This could seriously ruin your files if you use it
%incorrectly

% Path definitions
%h5Path = '/Users/jmckenzi/Desktop/LUT/';
h5Path = '/Volumes/Data/Data/Ovarian/immuno_Tcells/39 samples/';
xlPath = '/Users/jmckenzi/Desktop/LUT/LUT.xlsx';

% Read in the xl file
[~,~,xl] = xlsread(xlPath);

% Which h5 files are in the folder?
files = fileFinderAll(h5Path,'h5');
files = files(2:37,:);
numF = size(files,1);

% Loop
for n = 1:numF
    
    % Find a direct filename match in the spreadsheet
    fx = strcmp(xl(:,1),files{n,2}(1:end-3));
    if sum(fx) ~= 1
        disp('Error');
        continue;
    end
    
    % Read the annoatation parts from the h5 file...
    info = h5info([files{n,1} filesep files{n,2}],'/tissue_id')
    
    % Read in all the attributes, i.e. annotations...
    numA = size(info.Attributes,1) - 1;
    atts = cell(numA,1);
    for r = 1:numA
        atts{r,1} = lower(h5readatt([files{n,1} filesep files{n,2}],...
            '/tissue_id',int2str(r)));
    end
    
    % Now we need to do the specified replacement
    stroma = ~cellfun(@isempty,strfind(atts,'stroma'));
    if sum(stroma) == 1
        atts(stroma) = xl(fx,15);
    end
    
    % Now tumour
    tumour = ~cellfun(@isempty,strfind(atts,'tumour')) | ...
        ~cellfun(@isempty,strfind(atts,'tumor'));
    if sum(tumour) == 1
        atts(tumour) = xl(fx,16);
    end
    
    atts
    
    % Now just write them out to the original file
    for r = 1:numA
        h5writeatt([files{n,1} filesep files{n,2}],'/tissue_id',...
            int2str(r),atts{r,1});
    end
    
    
end




end


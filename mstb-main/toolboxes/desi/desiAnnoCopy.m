function desiAnnoCopy
%desiAnnoCopy - copy annotations from one mat file to another. File names
%must be the same for this to work

annoPath = '/Volumes/Data/Data/PBC Liver/MAT files/';
blankPath = '/Users/jmckenzi/DB/PBC/';

% Get the file list for both locations
annoList = fileFinderAll(annoPath,'mat');
blankList = fileFinderAll(blankPath,'mat');

% Loop through the blank files looking for matches in the anno list
for n = 2:size(blankList,1)
    
    % Find matching file names
    fx = strcmp(annoList(:,2),blankList{n,2});
    
    % If there is a single match, then we can print them
    if sum(fx) == 1
        aF = [annoList{fx,1} filesep annoList{fx,2}];
    else
        continue;
    end
    
    % Now open aF and get the annotations and optical image parts
    tmp = open(aF);
    
    % Open the blank file into which we will copy the annotations etc
    blf = open([blankList{n,1} filesep blankList{n,2}]);
    blf = blf.dpn;
    
    % Copy in the parts...
    blf.opt = tmp.dpn.opt;
    blf.anno = tmp.dpn.anno;
    dpn = blf;
    
    % Save again...
    save([blankList{n,1} filesep blankList{n,2}],'dpn');
    
    disp([int2str(n) '/' int2str(size(blankList,1))]);
end


end


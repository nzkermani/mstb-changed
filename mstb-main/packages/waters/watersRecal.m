function [ orig ] = watersRecal(raw,replace)
%watersRecal - for a given centroid mode file, replace the cal line
%function either with the 000 line, or with the original text in order to
%restore the file as it was

% From the 'raw' file path, determine the name of the header file
head = [raw filesep '_HEADER.TXT'];

% What are we looking for?
calLine = '$$ Cal Function 1:';

% Cell of arbitrary length
txt = cell(100,1);
n = 0;

% Create empty lineID variable
lineID = 0;

% Read in the file
fid = fopen(head,'r');
while ~feof(fid)    
    n = n + 1;
    txt{n,1} = fgetl(fid);
    
    % Identify the correct line...
    if length(txt{n,1}) >= length(calLine)
        if strcmp(txt{n,1}(1:length(calLine)),calLine)
            lineID = n;
        end
    end
       
end
fclose(fid);

% If the lineID remains at 0, then we have not found the appropriate line
% in the folder. As such, we just quit...
if lineID == 0
    orig = '';
    return;
end

% Trim the list a little
txt = txt(1:n,:);

% Now we have read in the file, we can decide if we want to replace the
% line with the 000 or the original...
if length(replace) == 0 %#ok<ISMT>
    % Then we need to dump in the zeros as suggested by Keith
    newLine = [calLine ' 0.0,1.0,0.0,0.0,0.0,0.0,T1'];
    orig = txt{lineID,1};
    
    % Replace this newLine in the txt cell array
    txt{lineID,1} = newLine;

    % Set the flag!
    flag = 'orig2new';
    
else
    % Then restore what it was before
    flag = 'restore';
end


% Replace or restore the original file
switch flag
    
    case 'orig2new'
        
        % Move the old file
        copyfile(head,[raw filesep 'ORIG_HEADER.TXT'],'f');
        
        % Now write the new one
        fid = fopen(head,'w');        
        for n = 1:size(txt,1)
            fprintf(fid,'%s\n',txt{n,1});
        end
        fclose(fid);
        
    case 'restore'
        
        % Restore from the old file name
        copyfile([raw filesep 'ORIG_HEADER.TXT'],head,'f');
        orig = [];
        
end
   
end


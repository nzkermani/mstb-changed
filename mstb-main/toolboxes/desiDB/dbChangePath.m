function dbChangePath(src,event,fig)
%dbChangePath - alter the main location of the DB

% See if there is already a path that we could use from the table
try
    tab = get(fig.tab,'Data');
    defP = tab{1,2};
catch
    defP = pwd;
end

defP = uigetdir(defP);
if isnumeric(defP)
    return
end
defP = [defP filesep];

% Now run the function to determine the files that will go in the database
[allF] = dbFileFind(defP);

% Get the file information for the table
[tabDat] = dbFileInfo(allF);

% Set the table data
dbTableUpdate(fig,tabDat);
set(fig.fig,'Name',defP);


end


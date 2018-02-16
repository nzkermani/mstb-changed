function pdCallbacks(fig,defP)
%pdCallbacks - add the callbacks for the various parts of pd

% File selection...
set(fig.fSelect,'Callback',{@pdFileSelect,fig,defP});

% Switch between pos / neg
set(fig.polarity,'Callback',{@pdPolaritySwitch,fig});

% Extract ions
set(fig.extract,'Callback',{@pdIonExtract,fig});

% Draw image - can include all visualise controls with same function
set(fig.ionDisplay, 'Callback',{@pdDisplay,fig});
set(fig.minInt,     'Callback',{@pdDisplay,fig});
set(fig.maxInt,     'Callback',{@pdDisplay,fig});
set(fig.doLog,      'Callback',{@pdDisplay,fig});

% Determine the number of segments, or the location to be more accurate
set(fig.doSegment,'Callback',{@pdSegmentDetermine,fig});

% Delete segment button
set(fig.segDelete,'Callback',{@pdSegmentDelete,fig});

% Add segment button
set(fig.segAdd,'Callback',{@pdSegementAdd,fig});

% Process button
set(fig.process,'Callback',{@pdProcess,fig});

end


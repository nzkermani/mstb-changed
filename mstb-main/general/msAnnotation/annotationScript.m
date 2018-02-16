%% Annotation Script
% Let's you do annotation and identification of specific trends within a
% dataset. There are a lot of individual sections that need to be run.
%
% You will need to have the following as inputs:
% mz        - vector of m/z values
% sp        - matrix of spectral intensities
% histID    - cell vector of histological identities

% First annotate the variables
[lm,ass] = annotateMZ(mz,'ppm',10);

% Now determine the lipid class of variables that differ significantly
% between histological groups. list is a cell array containing all of the
% information about the annotated lipid, e.g. index|m/z|p|q|grp|lipid|class
[list] = annotateLipidClass(mz,sp,histID,ass,lm);

% This decides on the annotation factors and returns a useful list showing
% variables and their characterisitcs. The annoStr structure contains
% tables showing the number of features that are significant for each group
[annoInfo,annoStr] = annotateByGroup(list)

% So then finally we need some way to plot figures to show which groups etc
% are significant
annotateFigure(list,annoInfo,annoStr,type);

%% Example script for annotation of spectral data.
% For instructions and a guide, please refer to the documentation available
% here: 

% Initial annotation
[lm,ass] = annotateMZ(op.cmz,...
    'Polarity','negative',...
    'Adduct',{'M-H','M+Cl'},...
    'Tolerance',10,...
    'Database','luisa');

%% Exporting the annotations
qValueThreshold = 2;
annotateOP2(op.cmz,full(op.XPeaks),op.histID,ass,lm,...
    'Test2.txt',...
    qValueThreshold);

%% Annotation figure
annotateClassInfo(lm,ass);

%% Significant lipids
[list] = annotateLipidClass(mz,sp,...
    histID,ass,lm);

%% Annotate by group
[factDat,annoInfo] = annotateByGroup(list);

%% Annotation figures
annotateFigure(list,factDat,...
    annoInfo,'pi');

annotateFigure(list,factDat,...
    annoInfo,32);

annotateFigure(list,factDat,...
    annoInfo,'class');
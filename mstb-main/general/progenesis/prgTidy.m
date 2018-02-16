function [ x ] = prgTidy(op)
%prgTidy - sort out the progenesis import function because there is a lot
%of extra rubbish that gets in the way...


%% VARIABLE INFORMATION
% Convert mz/rt to numbers - the rest is not important
tmp = strcmp(op.varH,'m/z');
mz = str2double(op.var(:,tmp));
tmp = strcmp(op.varH,'Retention Time');
rt = str2double(op.var(:,tmp));

% Put in structure
tmp = strcmp(op.varH,'Feature Name');
x.v.rtmz = op.var(:,tmp);
x.v.rt = rt;
x.v.mz = mz;
x.v.all = op.var;
x.v.allH = op.varH;

%% Y DATA / METADATA
% Extract some of the more important values...
tmp = strcmp(op.metaH,'Run Order');
x.y.meta.run = str2double(op.meta(:,tmp));
tmp = strcmp(op.metaH,'Plate');
x.y.meta.plate = str2double(op.meta(:,tmp));
tmp = strcmp(op.metaH,'Batch');
x.y.meta.batch = str2double(op.meta(:,tmp));

tmp = strcmp(op.metaH,'Patient ID');
x.y.meta.patientID = op.meta(:,tmp);
tmp = strcmp(op.metaH,'Tissue');
x.y.meta.tissueID = op.meta(:,tmp);
x.y.all = op.meta;
x.y.allH = op.metaH;

%% X DATA
x.x = op.data;

end


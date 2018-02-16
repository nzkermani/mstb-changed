 %% example to run the workflow
 % pathName: directory of the imzMl files
 % sample: name of the imzML sample without any extension (Please make sure
 % that h5 file present in the pathName directory with the same sample name(sample.h5).
 % this is necessary for the test for spatial randomness)
 PathName = '/Users/jmckenzi/DB/ColorectalDummy/A30/';
 sample = 'A30centreS475um';
 [SPnew, MZnew, bin, SizeMZ]=workflowDESI(PathName,sample);
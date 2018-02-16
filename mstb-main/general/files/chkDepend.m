function [ fList ] = chkDepend(mFile)
%chkDepend - find dependent files 

clear functions

[fList,pList] = matlab.codetools.requiredFilesAndProducts(mFile);

disp('%%%%%%%%%%%%% DEPENDENCIES %%%%%%%%%%%%%');
disp(fList')
disp([char(8) '%%%%%%%%%%%%% DEPENDENCIES %%%%%%%%%%%%%' char(10)]);

end


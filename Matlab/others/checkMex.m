function checkMex()
% compile cpp files if neccesary. You need to set up your compiler
% configuration with "mex -setup" beforehand.
if exist(fullfile('.','kernel',['poly2d.' mexext]),'file')
    disp(['Compiled file poly2d.' mexext ' exist!']);
else
    curPath=cd;
    disp('Compiling ply2d.cpp ... ...')
    cd(fullfile('.','kernel'));
    mex('poly2d.cpp');
    cd(curPath);
end

% add code for compling other cpp files below ... ...

%提交INP文件,测试成功
%输入：Path，INP文件所在的绝对位置
%输入：Inpfile,Inp文件名，不要包含后缀名
%输入：cpus,调用cpu个数，需要引号括起来
function [outputArgs] = submitINP(Path,InpFile,cpus)
inputFile=[' abaqus job=',InpFile,' cpus=',cpus];                           %无用户子程序
t0=tic;                                                                     %开始计时
MatlabPath=pwd();                                                           %记下当前matlab目录
cd(Path);                                                                   %进入Abaqus目录,即inp文件所在目录
[outputArgs] =system(inputFile);                                            %通过系统调用，运行ABAQUS,提交计算文件
pause(5)
cd(MatlabPath);                                                             %返回Matlab工作目录
if (exist([Path,'\',InpFile,'.lck'],'file')==2)                             %若提交成功，则检测计算时间
    while exist([Path,'\',InpFile,'.lck'],'file')==2
        t=toc(t0);
        h=fix(t/3600);
        m=fix(mod(t,3600)/60);
        sec=fix(mod(mod(t,3600),60));
        pause(1)
        fprintf(' ----------ABAQUS calculating----------\n 花费时间 %d:%d:%d\n',...
            h,m,sec);
    end
    fprintf('----------ABAQUS 计算完成----------\n 花费时间 %d:%d:%d\n', ...
        h,m,sec);
else
    fprintf('\n INP文件提交失败 \n')
end
end

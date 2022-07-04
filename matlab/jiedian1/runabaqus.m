%向abaqus提交计算文件,测试成功
%输入：Path，inp计算文件所在的绝对路径。
%输入：UserFile：用户Fortran子程序，如果有子程序就给子程序的文件名，没有的话就不要这个参数，
% 把runabaqus的一行代码
% inputFile=[‘abaqus job=’,InpFile,’ user=’,userFile,’ cpus=’,cpus]，
% 改为inputFile=[‘abaqus job=’,InpFile,’ cpus=’,cpus]。
%输入：InpFile，inp计算文件的文件名
%输入：cpus,指定ABAQUS的求解器CPU数量
function [outputArgs] = runabaqus(Path,InpFile,cpus)
inputFile=['cmd/c abaqus job=',InpFile,' cpus=',cpus,' int ask=off'];       %无用户子程序
t0=tic;                                                                     %开始计时
MatlabPath=pwd();                                                           %记下当前matlab目录
cd(Path);                                                                   %进入Abaqus目录,即inp文件所在目录
[outputArgs] =system(inputFile);                                            %通过系统调用，运行ABAQUS,提交计算文件
pause(5)
cd(MatlabPath); %返回Matlab工作目录
% if (exist([Path,'\',InpFile,'.lck'],'file')==2)    
%     %若提交成功，则检测计算时间
%     while exist([Path,'\',InpFile,'.lck'],'file')==2
%         t=toc(t0);
%         h=fix(t/3600);
%         m=fix(mod(t,3600)/60);
%         sec=fix(mod(mod(t,3600),60));
%         pause(1)
%         fprintf(' ----------ABAQUS calculating----------\n 花费时间 %d:%d:%d\n',...
%         h,m,sec);
%     end
%     fprintf('----------ABAQUS 计算完成----------\n 花费时间 %d:%d:%d\n', ...
%     h,m,sec);
% else
%     %若提交计算出错,则输出错误信息
%     fprintf('\n INP文件提交失败 \n')
% end
end

%运行python脚本，读取odb数据库的历史输出
function getHistoryOutput(Path,OdbFile,step,req)
% %specified part instance,node set,output
%参数写入Abaqus工作目录
ReqFile=[Path,'\req.txt'];
fid=fopen(ReqFile,'wt');
fprintf(fid,'%s,%s,%s,%s',Path,OdbFile,step,req);
fclose(fid);
%写入当前目录
ReqFile='req.txt';
fid=fopen(ReqFile,'wt');
fprintf(fid,'%s,%s,%s,%s',Path,OdbFile,step,req);
fclose(fid);
%execute python file
system('abaqus cae noGUI=odbHistoryOutput.py');%通过dos命令调用python脚本
showlogfile('pylog.txt');%显示python运行的输出信息
end
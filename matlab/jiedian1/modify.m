function modify(Path,InpFile,NewData)
%modify.m
%libing403,2017-05-09
ReqFile=[Path,'\modify.txt'];
fid=fopen(ReqFile,'wt');
%把需要修改的inp文件和新数据文件的文件名、路径写入modify.txt
fprintf(fid,'%s,%s.inp,%s.txt',Path,InpFile,NewData);
fclose(fid);
system('abaqus cae noGUI=inpmodify.py');%调用python脚本修改inp数据
end

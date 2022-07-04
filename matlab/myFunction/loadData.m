%将其他目录下其他格式的文件转换成mat格式并加载到工作区
%可以加载的文件格式：txt,elsx,mat,dat
function [returnvalue]=loadData(filePath,fileName,fileFormat)
MatlabPath=pwd();%记下当前matlab目录
cd(filePath);
if (exist([filePath,'\',fileName,'.',fileFormat],'file')==2)
    fprintf("存在文件，正在加载......\n");
    if string(fileFormat)=="txt"||string(fileFormat)=="mat"||string(fileFormat)=="xlsx"||string(fileFormat)=="dat"
        returnvalue=readmatrix([fileName,'.',fileFormat]);
        cd(MatlabPath); %返回Matlab工作目录 
        fprintf("加载完成！！！\n");
    else
        fprintf("暂不支持该文件类型！！！\n");
    end
else 
    cd(filePath);
    fprintf("没有找到文件,请确认文件是否存在！");
end
end
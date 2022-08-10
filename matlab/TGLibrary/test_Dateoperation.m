%数据操作命令练习
% fopen:打开文件，返回文件标识符，保留0 1 2,输入参数r,w,a,a用于追加到现有文档的末尾
% fprintf：写入文件数据
% fclose
% ftell
% fgetl
% ftell:返回光标当前位置
% fseek
% fwrite；写入二进制文件数据
%%
clc;clear;close all
fid = fopen('E:\code\matlab\TGLibrary\myData\badpoem.txt');%打开文件
ftell(fid)%查询当前位置
tline1 = fgetl(fid) ; %读取一行的数据
ftell(fid)%查询当前位置，返回
tline2 = fgetl(fid);
ftell(fid)%查询当前位置
tline3 = fgetl(fid);
ftell(fid)%查询当前位置
fseek(fid,20,'bof');%移动文件到指定的位置，bof指定文件的开头
fgetl(fid);%查询当前位置
fclose(fid); %关闭文件
%%
clc;clear;close all
x = 0:0.1:1;
y = [x; exp(x)];
fileID = fopen('E:\code\matlab\TGLibrary\myData\exptable.txt','w');
fprintf(fileID, 'Exponential Function\n\n');
fprintf(fileID,'%f %f\n',y);
fclose(fileID);
type E:\code\matlab\TGLibrary\myData\exptable.txt
%%
clc;clear;close all
fileID = fopen('nine.bin','w');
fwrite(fileID,[1:9]);
fclose(fileID);





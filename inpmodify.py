#! /user/bin/python
#- -coding: UTF-8-*-
#inpmodify.py
# 本脚本的功能是修改ABAQUS inp文件
#libing403，2017-5-7
import time
#读取inp文件的路径及文件名内容
f=open('modify.txt','r')
req=f.readline()
f.close()
req=req.split(',')
InpFile=req[0]+'/'+req[1]
NewData=req[0]+'/'+req[2]
#读取inp文件内容
fid=open(InpFile,"r")
lines=fid.readlines()
fid.close()
#找出原来的数据行
startstr="*Permeability, specific=1.\n"
startIndex=lines.index(startstr)+1
#28个行数据需要替换
endIndex=startIndex+28
#读入新数据
fid=open(NewData,"r")
newInp=fid.readlines()
fid.close()
print("%s"%newInp)
#替换原来的数据行
i=0
for Index in range(startIndex,endIndex):
    lines[Index]=newInp[i]
    i=i+1
#写入新数据，替换原来的数据
fid=open(InpFile,"w")
fid.writelines(lines)
fid.close()
#写入操作日志
meg="inpmodify message:\n inp file modify successfully\n"
fid=open("pylog.txt","a")
fid.write('%s\n'%meg)
fid.close()
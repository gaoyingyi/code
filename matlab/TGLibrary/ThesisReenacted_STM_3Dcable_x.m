%[1]王勇,高冀峰.桁架结构几何非线性问题的力密度法[J].应用力学学报,2011,28(05):509-513+556.
%文献复演
clc;clear;close all
syms xi xj yi yj zi zj E A L0
www=[xi yi zi xj yj zj].';
L=sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);%变形后的长度L
N=E*A*(L-L0)/L0;%原长L0，假设已知，是一个常量
fxi=N*(xi-xj)/L;
fxj=-N*(xi-xj)/L;
fyi=N*(yi-yj)/L;
fyj=-N*(yi-yj)/L;
fzi=N*(zi-zj)/L;
fzj=-N*(zi-zj)/L;
fai=[fxi,fyi,fzi,fxj,fyj,fzj].';
K=jacobian(fai,www);
K=subs(K,str2sym('sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2)'),str2sym('L'));
K=subs(K,str2sym(' xi - xj'),str2sym('u'));
K=subs(K,str2sym(' yi - yj'),str2sym('v'));
K=subs(K,str2sym(' zi - zj'),str2sym('w'));
K=simplify(K);
K=expand(K)
clearvars -except K;
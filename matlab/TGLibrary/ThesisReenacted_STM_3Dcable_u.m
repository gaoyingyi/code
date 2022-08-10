%[1]王勇,高冀峰.桁架结构几何非线性问题的力密度法[J].应用力学学报,2011,28(05):509-513+556.
%文献已经在其他文件中复演。
%论文存在一定问题：
% 原论文推导的刚度矩阵是关于坐标的函数，而不是位移；
% 原论文令变形后的长度是L，而原长L0假设已知，是一个常量；
%本文件推导刚度矩阵关于位移的函数；
%原长L0是关于节点坐标的常量
%变形后的长度L是位移的函数
%本文件推导的结果与以下论文基于工程应变泰勒展开推导的结果基本一致，P62
%[2]刘望. 星载大型索网天线结构非线性动力学建模及形面误差特性研究[D].国防科学技术大学,2013.
%所不同的是，刘望，将L泰勒展开，本文件推导的结果直接代入，
%推导方法也不同。最后结果只在形式上有差异。刘望文中将非线性项细分，这里不再展开；
clc;clear;close all
syms xi xj yi yj zi zj ui uj vi vj wi wj N0 E A L0
u=[ui vi wi uj vj wj].';%节点位移向量
L=sqrt((xj+uj-xi-ui)^2+(yj+vj-yi-vi)^2+(zj+wj-zi-wi)^2);
N=N0+E*A*(L-L0)/L0;
fxi=-N*(xj+uj-xi-ui)/L;
fxj=N*(xj+uj-xi-ui)/L;
fyi=-N*(yj+vj-yi-vi)/L;
fyj=N*(yj+vj-yi-vi)/L;
fzi=-N*(zj+wj-zi-wi)/L;
fzj=N*(zj+wj-zi-wi)/L;
fai=[fxi,fyi,fzi,fxj,fyj,fzj].';
K=jacobian(fai,u);

K=subs(K,str2sym('sqrt((xj+uj-xi-ui)^2+(yj+vj-yi-vi)^2+(zj+wj-zi-wi)^2)'),str2sym('L'));
K=subs(K,str2sym('sqrt((xj-xi)^2+(yj-yi)^2+(zj-zi)^2)'),str2sym('L0'));
%K=subs(K,str2sym('N0+E*A*(L-L0)/L0'),str2sym('N'));
K=simplify(K);
K=subs(K,str2sym('xj+uj-xi-ui'),str2sym('u*L'));
K=subs(K,str2sym('yj+vj-yi-vi'),str2sym('v*L'));
K=subs(K,str2sym('zj+wj-zi-wi'),str2sym('w*L'));
K=simplify(K);
K=expand(K)
save('E:\code\matlab\TGLibrary\myData\Kcable.mat',"K");
clearvars -except K;
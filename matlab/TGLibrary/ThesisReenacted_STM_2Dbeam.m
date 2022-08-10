%周水兴,陈山林.MatLab在有限元刚度矩阵推导中的应用[J].重庆交通学院学报,2007(02):29-31.
%二维梁刚度矩阵推导
clc;clear;close all;
syms x L y E h W A I
Nu=[1-x/L, 0, 0, x/L,0, 0];%定义形函数矩阵Nu
Nv=[0, 1-3*x^2/L^2+2*x^3/L^3,x-2*x^2/L+x^3/L^2,0,3*x^2/L^2-2*x^3/L^3,-x^2/L+x^3/L^2];
B=[diff(Nu,'x',1);-y*diff(Nv,'x',2)];
D=[E,0;0,E;];
K=B.'*D*B;
K=int(K,'x','0','L');
K=int(K,'y',str2sym('-h/2'),str2sym('h/2'));
K=int(K,'K',str2sym('-W/2'),str2sym('W/2'));
K=subs(K,str2sym('h^3*W'),str2sym('12*I'));
K=subs(K,str2sym('h*W'),str2sym('A'));
K=simplify(K)
clearvars -except K




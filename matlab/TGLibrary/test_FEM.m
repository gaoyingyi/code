% clc;clear;close all
% E=210000;%弹性模量
% mu=0.3;%泊松比
% N0=10;%初轴向力
% r=5;%横截面半径
% rho=1.2e-5;
% nodeCoordinateMatrix=[100,-300,500;1000,2000,500;];%全局坐标系下节点坐标矩阵
% nodeDisplacementMatrix=[0.11,0.33,0.5,0,0,0.2,0,0,0.8,0,0.9,0].';%局部坐标系下节点位移
% Kbeam=ESM_B31(E,mu,r,N0,nodeCoordinateMatrix,nodeDisplacementMatrix);
% Mbeam=EMM_B31(E,mu,rho,r,nodeCoordinateMatrix);
%% newton-Raphson求解非线性平衡方程的静力位移问题
clc;clear;close all 
syms x y z fx1 fy1 fz1
E=210000;%弹性模量
mu=0.3;%泊松比
N0=0.0001;%初轴向力
r=5;%横截面半径
rho=1.2e-5;
nodeCoordinateMatrix=[501,501,501;1020,1020,1000;];%全局坐标系下节点坐标矩阵
nodeDisplacementMatrix=[0,0,0,x,y,z].';%全局坐标系下节点位移
u0=[0 0 0].';
FTT=[fx1,fy1,fz1,1,1,1].';
Ku=ESM_cable(E,mu,r,N0,nodeCoordinateMatrix,nodeDisplacementMatrix);
K=Ku(4:6,4:6);
u=nodeDisplacementMatrix(4:6,:);
F=FTT(4:6,:);
Fai=K*u-F;
solv=NR_Method(Fai,K,u0)
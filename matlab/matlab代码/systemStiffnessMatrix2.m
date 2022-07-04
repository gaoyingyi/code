%加上非线性弹簧的系统模型测试
clc;clear;
%% 参数
T=0;%索的内拉力
F=0;%梁的压力
LCable=400;
Lbeam=400;

%% 局部坐标系下单元刚度矩阵
[kspring]=springElementStiffnessMatrix(1);%局部坐标系弹簧单元矩阵
[kcable,mcable]=cableElementStiffnessMatrix(T,LCable);%局部坐标系索单元矩阵
[kbeam,mbeam]=beamElementStiffnessMatrix(F,Lbeam);%局部坐标系梁单元矩阵
%% 整体坐标系下单元刚度矩阵
springCTF=CTF13(400,0,0,320,60,0);%弹簧的转换矩阵
kspring=springCTF'*kspring*springCTF;%整体坐标系下弹簧的刚度矩阵
kspring=[kspring(1:3,1:3),zeros(3),kspring(1:3,4:6),zeros(3);
         zeros(3),zeros(3),zeros(3),zeros(3);
         kspring(4:6,1:3),zeros(3),kspring(4:6,4:6),zeros(3);
         zeros(3),zeros(3),zeros(3),zeros(3);];%扩充成12×12方便组集

lambdaC=coordinateTransformation(320,60,0,0,300,0,0,0,0);
cableCTF=[lambdaC,zeros(3),zeros(3),zeros(3);
         zeros(3),lambdaC,zeros(3),zeros(3);
         zeros(3),zeros(3),lambdaC,zeros(3);
         zeros(3),zeros(3),zeros(3),lambdaC;];%索的转换矩阵
kcable=cableCTF'*kcable*cableCTF;%整体坐标系下索的刚度矩阵

lambdaB=coordinateTransformation(0,0,0,400,0,0,0,300,0);
beamCTF=[lambdaB,zeros(3),zeros(3),zeros(3);
         zeros(3),lambdaB,zeros(3),zeros(3);
         zeros(3),zeros(3),lambdaB,zeros(3);
         zeros(3),zeros(3),zeros(3),lambdaB;];%梁的转换矩阵
kbeam=beamCTF'*kbeam*beamCTF;%整体坐标系下梁的刚度矩阵
%% 扩充弹簧单元刚度矩阵
Z=zeros(24);
Z(7:18,7:18)=kspring;
kspring=Z;
clear Z;
%% 扩充索单元刚度矩阵
Z=zeros(24);
Z(13:24,13:24)=kcable;
kcable=Z;
clear Z;
%% 扩充梁单元
Z=zeros(24);
Z(1:12,1:12)=kbeam;
kbeam=Z;
clear Z;
%% 集成总刚
K=kspring+kbeam+kcable;
%% Z的负方向加10N的力验证模型的正确性
u=K(7:14,7:14)\[0,-10,0,0,0,0,0,0]';
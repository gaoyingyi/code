%系统模型
clc;clear;
%参数
T=0;%索的内拉力
LCable=500;
Lbeam=400;
F=0;%梁的压力
%% 处理弹簧单元
Kspring=500; %定义非线性刚度
ks=[1,0,0,-1,0,0;
   0,0,0,0,0,0;
   0,0,0,0,0,0;
  -1,0,0,1,0,0;
   0,0,0,0,0,0;
   0,0,0,0,0,0;];
Kspring=Kspring.*ks;%扩充三个方向后的弹簧单元刚度矩阵
lambdaA=coordinateTransformation(320,60,0,400,0,0,320,0,0);
springCTF=[lambdaA,zeros(3);
          zeros(3),lambdaA];
Kspring=springCTF'*Kspring*springCTF;%整体坐标系下弹簧的刚度矩阵
u=inv(Kspring(4:6,4:6))*[10,0,0]'%经过验证，与abaqus计算模型基本一致。
%% 局部坐标系下单元刚度矩阵
[kcable,mcable]=cableElementStiffnessMatrix(T,LCable);
[kbeam,mbeam]=beamElementStiffnessMatrix(F,Lbeam);
%% 整体坐标系下单元刚度矩阵
lambdaC=coordinateTransformation(0,300,0,320,60,0,0,0,0);
cableCTF=[lambdaC,zeros(3);
          zeros(3),lambdaC];%索的转换矩阵
kcable=cableCTF'*kcable*cableCTF;%整体坐标系下索的刚度矩阵
lambdaB=coordinateTransformation(0,0,0,400,0,0,0,300,0);
beamCTF=[lambdaB,zeros(3),zeros(3),zeros(3);
         zeros(3),lambdaB,zeros(3),zeros(3);
         zeros(3),zeros(3),lambdaB,zeros(3);
         zeros(3),zeros(3),zeros(3),lambdaB;];%梁的转换矩阵
kbeam=beamCTF'*kbeam*beamCTF;%整体坐标系下梁的刚度矩阵
%% 扩充索单元刚度矩阵
Z=zeros(15);
Z(7:12,7:12)=kcable;
kcable=Z;
%% 扩充梁单元
Z=zeros(15);
Z(1:6,1:6)=kbeam(1:6,1:6);
Z(1:6,10:15)=kbeam(1:6,7:12);
Z(10:15,1:6)=kbeam(7:12,1:6);
Z(10:15,10:15)=kbeam(7:12,7:12);
kbeam=Z;
clear Z;
%% 集成总刚
K=kbeam+kcable;
%% Z的负方向加10N的力验证模型的正确性
u=inv(K(10:15,10:15))*[10,-10,10,0,0,0]'%经过验证，与abaqus计算模型基本一致。
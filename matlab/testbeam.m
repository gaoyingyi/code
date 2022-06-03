clc;clear;
%% 梁的结果测试，x方向加100N外力，与手算结果一致，但是内力无影响。刚度矩阵哪=那里的问题，就是无影响。
[kbeam,mbeam]=beamElementStiffnessMatrix(100,1000);
syms x2 y2 z2 theta1 theta2 theta3 f1 f2 f3 t1 t2 t3
ue=[0,0,0,0,0,0,x2,y2,z2,theta1,theta2,theta3]';
F=[f1,f2,f3,t1,t2,t3,100,0,0,0,0,0]';
u=inv(kbeam(7:12,7:12))*[100,0,0,0,0,0]';%计算的位移与未加应力的位移相同

%% 索的位移测试，不加应力计算结果与手算相同，加上应力位移有变化。
% [kbeam,mbeam]=cableElementStiffnessMatrix(100,1000);
% syms x2 y2 z2 theta1 theta2 theta3 f1 f2 f3 t1 t2 t3
% ue=[0,0,0,0,0,0,x2,y2,z2,theta1,theta2,theta3]';
% F=[f1,f2,f3,t1,t2,t3,100,0,0,0,0,0]';
% u=inv(kbeam(4:6,4:6))*[100,0,0]';




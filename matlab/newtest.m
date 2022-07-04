clc;clear;close all;
%%参数定义
global sumsystemDof sumBoundarynodeDof EBeam muBeam rhoBeam ECable muCable rhoCable massDampingCoefficient stiffnessDampingCoefficient;
sumsystemDof=18;                                          %系统总的自由度数目
sumBoundarynodeDof=9;                                     %边界节点自由度数目
EBeam=350000;                                                    %梁弹性模量
muBeam=0.3;                                                        %梁泊松比
rhoBeam=1.5e-9;                                                      %梁密度
ECable=20000;                                                    %梁弹性模量
muCable=0.3;                                                       %梁泊松比
rhoCable=1.5e-9;                                                     %梁密度
massDampingCoefficient=0;                                      %质量阻尼系数
stiffnessDampingCoefficient=0;                                 %刚度阻尼系数
%%
T=1;
deltaT=0.001;
[U,UD,UDD,sumQ]=newMark(T,deltaT); 

% 绘图

matlabU = U(2,:);
matlabUD = UD(2,:);
matlabUDD = UDD(2,:);

load("BBB.txt")
ttt=BBB(:,1)';
abaqusU=BBB(:,3)';
abaqusUD=BBB(:,4)';
abaqusUDD=BBB(:,2)';

tt = 0:deltaT:T;%绘图时间步

figure(1) 
hold on
    plot(tt,matlabU,'r','linewidth',0.1)
    plot(ttt,abaqusU,'b','linewidth',0.1)
    title('位移/转角')
figure(2)
hold on
    plot(tt,matlabUD,'r','linewidth',0.1)
    plot(ttt,abaqusUD,'b','linewidth',0.1)
    title('速度')
figure(3)
hold on
    plot(tt,matlabUDD,'r','linewidth',0.1)
    plot(ttt,abaqusUDD,'b','linewidth',0.1)
    title('加速度')

% for i=1:9
%     dn = U(i,:);    
%     dV = UD(i,:);    
%     dA = UDD(i,:);  
%     
%     figure(i)
%     subplot(3,1,1);
%     plot(tt,dn,'r','linewidth',0.1)
%     title('位移/转角')
%     subplot(3,1,2);
%     plot(tt,dV,'b','linewidth',0.1)
%     title('速度')
%     subplot(3,1,3);
%     plot(tt,dA,'g','linewidth',0.1)
%     title('加速度')
% end
% figure(10);
% Q=sumQ(2,:);%提取节点2 y方向上的集中力，即施加的外载荷，验证外载荷是否正确
% plot(tt,Q,'k','linewidth',0.1);title('外载荷'); 
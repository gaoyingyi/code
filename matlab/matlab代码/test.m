%% 测试力密度法，这是以前的测试，与有限元无关
% clc;clear;
% P=[0,0,0];
% C=[1,1,1,1,0,0,0,0]';
% Cf=[-eye(4);1,-1,0,0;0,-1,0,0;0,0,1,-1;-1,0,0,1];
% q=[2 2 2 2 1 2 3 4]';
% Q=diag(q);
% fixedNodeCoordinates=[100,0,0;0,100,0;-100,0,0;0,-100,0;];
% [x,y,z]=ForceDensityMethod(Q,C,Cf,P,fixedNodeCoordinates);
%% 测试newMark法,计算书上的例子，测试成功
% clc;clear;close all;
% NFreeDof=3;%除去位移边界后系统总的自由度数
% deltaT=0.363;%时间步长
% T=3.63; %总时间
% F=@(t) force(t,3);%调用外载荷力向量的函数句柄
% Ma=[1,0,0;
%    0,3,0;
%    0,0,1;];
% Ka=[2,-1,0;
%   -1,4,-2;
%    0,-2,2;];
% Ca=zeros(3);
% [U,UD,UDD,sumQ]=newMark(Ma,Ka,Ca,F,deltaT,T);
%% 测试网上找的悬臂梁,有初始速度的自由振动，测试成功
% 与网上找的代码跑出来的结果完全一样，证明了newmark函数的正确性
% 要修改newmark函数中的模块把相关内容去掉注释化  
% clc;clear;close all;
% rho =1;A = 1;E = 1;I = 1;L = 1;c1 = 0;c2 = 0;Ne = 100;
% deltaT=0.001;%时间步长
% T=20; %总时间
% [Ma,Ka] = Beam3(rho,A,E,I,L/Ne,Ne+1);
% Ca = c1*Ma+c2*Ka;
% F = @(t) 0;  
% [U,UD,UDD]=newMark(Ma,Ka,Ca,F,deltaT,T);
% dn = U(Ne-1,:);     % newmark中点位移
% dV = UD(Ne-1,:);    % newmark中点位移
% dA = UDD(Ne-1,:);   % newmark中点位移
% dt = 1e-3;T = 20;tt = 0:dt:T;
% hold on
% plot(tt,dn,'r','linewidth',2)
% plot(tt,dV,'b','linewidth',2)
% plot(tt,dA,'g','linewidth',2)
% title('Response of Mid-Point')
% xlabel('t')
% ylabel('Displacement')
% legend('Newmark')
%% 测试网上找的悬臂梁加中点集中力,测试成功
% 与网上找的代码跑出来的结果完全一样，证明了newmark函数的正确性
% 在99自由度上加了100*sin(2*pi*5*t)的集中力，要修改newmark函数中的模块把相关内容去掉注释化  
% clc;clear;close all;
% rho =1;A = 1;E = 1;I = 1;L = 1;c1 = 0.1;c2 = 0.1;Ne = 100;n=2;
% deltaT=0.001;%时间步长
% T=10; %总时间
% [Ma,Ka] = Beam3(rho,A,E,I,L/Ne,Ne+1);
% Ca = c1*Ma+c2*Ka;
% F = @(t) force(t,200);  
% [U,UD,UDD]=newMark(Ma,Ka,Ca,F,deltaT,T);
% dn = U(Ne-1,:);     % newmark中点位移
% dV = UD(Ne-1,:);    % newmark中点位移
% dA = UDD(Ne-1,:);   % newmark中点位移
% xx = 1/Ne:1/Ne:1;
% subplot(311)
% plot(0:deltaT:T,U(Ne-1,:),'linewidth',2)
% xlabel('t')
% title('Displacement of mid-point')
% subplot(312)
% plot(0:deltaT:T,U(2*Ne-1,:),'linewidth',2)
% xlabel('t')
% title('Displacement of end-point')
% subplot(313)
% plot([0,xx],[0;U(1:2:2*Ne-1,end)],'linewidth',2)
% xlabel('x')
% title('Mode Shape at the end time')
%% 测试网上的例子，每节点两个自由度，跑出来与代码图像参数完全一样,测试成功
% 证明newmark函数的正确
% clc;clear;close all;
% c1=0;c2=0;sumNode=101;NodeDof=2;
% deltaT=0.001;%时间增量
% T=4;         %总时间
% E = 7e10;rho = 2700;L = 0.67;b = 0.062;h = 0.0025;A = b*h;I = b*h^3/12;nL = 100;l = L/nL;
% %单元矩阵
% Ke = E*I/l^3.*[12,6*l,-12,6*l;
%              6*l,4*l^2,-6*l,2*l^2;
%              -12,-6*l,12,-6*l;
%              6*l,2*l^2,-6*l,4*l^2;];
% Me = rho*A*l/420.*[156,22*l,54,-13*l;
%                  22*l,4*l^2,13*l,-3*l*l;
%                  54,13*l,156,-22*l;
%                  -13*l,-3*l*l,-22*l,4*l*l;];
% KK = zeros(2*(nL+1),2*(nL+1));
% MM = zeros(2*(nL+1),2*(nL+1));
%              
% for i = 1:nL
%     KK((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2) = KK((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2)+Ke;
%     MM((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2) = MM((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2)+Me;
% end
% 
% K0 = KK(3:end,3:end);
% M0 = MM(3:end,3:end);
% C0 = c1*M0+c2*K0;
% F = @(t) force(t,200);
% [U,UD,UDD,sumQ]=newMark(M0,K0,C0,F,deltaT,T);
% %绘图
% Qi=sumQ(199,:);
% 
% dn = U(199,:);
% dV = UD(199,:);
% dA = UDD(199,:);
% tt =linspace(0,T,T/deltaT+1);
% 
% figure(1)
% plot(tt,dn,'r','linewidth',0.1);title('位移');
% figure(2)
% plot(tt,dV,'b','linewidth',0.1);title('速度');
% figure(3)
% plot(tt,dA,'m','linewidth',0.1);title('加速度');
% figure(4)
% plot(tt,Qi,'k','linewidth',0.1);title('外载荷'); 
%% 测试自己构想的悬臂梁模型，100个单元,测试成功
%  恒力外载荷，加阻尼和不加阻尼与abaqus最后稳态时结果基本一致，abaqus动力隐式分析步自带数值控制参数alpha，所以就算不加阻尼也会慢慢结束振动
%  受迫振动的振动频率与物体的固有频率无关，只与外载频率有关
%  与abaqus算出来的结果有区别，在网上查到，如果施加正弦激励最终梁会以激励的频率振动，这与算出来的结果是相符的。
%  认为matlab算出来的结果是正确的，abaqus可能有问题
% clc;clear;close all;
% FN=0;                                    %梁内力
% L=400;                                   %梁长度
% sumNode = 101;                           %总节点数
% sumElement=100;                          %总单元数
% NodeDof=6;                               %节点自由度数
% NFreeDof=NodeDof*sumNode-6;              %施加位移边界条件后，系统总的自由度数目
% LElement=L/sumElement;                   %单元长度
% c1 = 0;c2 = 0.3;                         %阻尼参数
% delta=0.5;alpha=0.25;                    %积分精度和稳定性决定的参数
% deltaT=0.01;                             %时间步长
% T=30;                                    %总时间
% K=zeros(NodeDof*sumNode,NodeDof*sumNode);%设置整体刚度矩阵
% M=zeros(NodeDof*sumNode,NodeDof*sumNode);%设置整体质量矩阵
% %C=zeros(n*Ne,n*Ne);                    %整体刚度矩阵初始化
% [Ke,Me]=beamElementStiffnessMatrix(FN,LElement);%单元刚度矩阵、质量矩阵
% 
% for i=1:sumElement                       %以单元的个数循环，组建整体刚度矩阵和整体质量矩阵
%     K(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)=K(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)+Ke;
%     M(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)=M(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)+Me;
% end
% 
% C = c1*M+c2*K;
% 
% Ka=K(7:end,7:end);                       %除去位移边界的自由度之后的刚度矩阵
% Ma=K(7:end,7:end);                       %除去位移边界的自由度之后的质量矩阵
% Ca=C(7:end,7:end);                       %除去位移边界的自由度之后的阻尼矩阵
% 
% F = @(t) force(t,NFreeDof);
% [U,UD,UDD,sumQ]=newMark(Ma,Ka,Ca,F,deltaT,T);
% 
% dn = U(596,:);    
% dV = UD(596,:);    
% dA = UDD(596,:); 
% Q=sumQ(596,:);
% 
% tt =linspace(0,T,T/deltaT+1);
% figure(1);
% plot(tt,dn,'r','linewidth',2);
% figure(2);
% plot(tt,dV,'b','linewidth',2);
% figure(3);
% plot(tt,dA,'g','linewidth',2);
% figure(4);
% plot(tt,Q,'k','linewidth',0.1);title('外载荷'); 
%% 测试自己构想的悬臂梁模型，1个单元,测试成功
%  结果
clc;clear;close all;
FN=0;                                    %梁内力
L=400;                                   %梁长度
sumNode = 2;                             %总节点数
sumElement=1;                            %总单元数
NodeDof=6;                               %节点自由度数
NFreeDof=NodeDof*sumNode-6;              %施加位移边界条件后，系统总的自由度数目
LElement=L/sumElement;                   %单元长度
c1 = 0.628;c2 = 0.3;                     %阻尼参数
delta=0.5;alpha=0.25;                    %积分精度和稳定性决定的参数
deltaT=0.01;                             %时间步长
T=10;                                    %总时间
K=zeros(NodeDof*sumNode,NodeDof*sumNode);%设置整体刚度矩阵
M=zeros(NodeDof*sumNode,NodeDof*sumNode);%设置整体质量矩阵
%C=zeros(n*Ne,n*Ne);                    %整体刚度矩阵初始化
[Ke,Me]=beamElementStiffnessMatrix(FN,LElement);%单元刚度矩阵、质量矩阵

for i=1:sumElement                       %以单元的个数循环，组建整体刚度矩阵和整体质量矩阵
    K(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)=K(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)+Ke;
    M(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)=M(1+6*(i-1):1+6*(i-1)+11,1+6*(i-1):1+6*(i-1)+11)+Me;
end
C = c1*M+c2*K;

Ka=K(7:end,7:end);                       %除去位移边界的自由度之后的刚度矩阵
Ma=M(7:end,7:end);                       %除去位移边界的自由度之后的质量矩阵
Ca=C(7:end,7:end);                       %除去位移边界的自由度之后的阻尼矩阵

F = @(t) force(t,NFreeDof);
[U,UD,UDD,sumQ]=newMark(Ma,Ka,Ca,F,deltaT,T);
dn = U(2,:);    
dV = UD(2,:);    
dA = UDD(2,:); 
Q=sumQ(2,:);

tt =linspace(0,T,T/deltaT+1);
figure(1);
plot(tt,dn,'r','linewidth',2);
figure(2);
plot(tt,dV,'b','linewidth',2);
figure(3);
plot(tt,dA,'g','linewidth',2);
figure(4);
plot(tt,Q,'k','linewidth',0.1);title('外载荷'); 
%% 测试索梁弹簧模型
% % clc;clear;close all;
% %参数
% N=0;%梁内力
% TT=0;%索单元内力
% sumNode = 3;%总节点数
% everynodeDOF=6;%每节点自由度数
% sumBC_Displacement=9;%位移边界条件的个数
% Lbeam=400;Lcable=400;
% c1 =0.0625;c2=0.008;
% delta=0.5;alpha=0.25;
% deltaT=0.001;%时间步长
% T=100;     %总时间
% % 全局坐标系下的单元刚度
% K=zeros(everynodeDOF*sumNode,everynodeDOF*sumNode);%初始化全局刚度矩阵
% M=zeros(everynodeDOF*sumNode,everynodeDOF*sumNode);%初始化全局质量矩阵
% %C=zeros(everynodeDOF*sumNode,everynodeDOF*sumNode);%初始化全局阻尼矩阵
% 
% [Kebeam,Mebeam]=beamElementStiffnessMatrix(N,Lbeam);      %梁单元刚度质量矩阵
% [Kecable,Mecable]=cableElementStiffnessMatrix(TT,Lcable); %索单元刚度质量矩阵
% [Kespring,Mespring]=springElementStiffnessMatrix(1);      %弹簧单元刚度质量矩阵
% %全局坐标系下梁刚度矩阵，质量矩阵
% lambdaB=coordinateTransformation(0,0,0,400,0,0,0,300,0);
% beamCTF=[lambdaB,zeros(3),zeros(3),zeros(3);
%          zeros(3),lambdaB,zeros(3),zeros(3);
%          zeros(3),zeros(3),lambdaB,zeros(3);
%          zeros(3),zeros(3),zeros(3),lambdaB;];%梁的转换矩阵
% Kebeam=beamCTF'*Kebeam*beamCTF;%整体坐标系下梁的刚度矩阵
% Mebeam=beamCTF'*Mebeam*beamCTF;%整体坐标系下梁的刚度矩阵
% %全局坐标系下索刚度矩阵，质量矩阵
% lambdaC=coordinateTransformation(320,60,0,0,300,0,0,0,0);
% cableCTF=[lambdaC,zeros(3);
%           zeros(3),lambdaC;];%索的转换矩阵
% Kecable=cableCTF'*Kecable*cableCTF;%整体坐标系下索的刚度矩阵
% Mecable=cableCTF'*Mecable*cableCTF;%整体坐标系下索的刚度矩阵
% %全局坐标系下弹簧刚度矩阵，质量矩阵
% lambdaD=coordinateTransformation(400,0,0,320,60,0,0,180,0);%弹簧的转换矩阵
% springCTF=[lambdaD,zeros(3);
%            zeros(3),lambdaD;];%索的转换矩阵
% Kespring=springCTF'*Kespring*springCTF;%整体坐标系下索的刚度矩阵
% Mespring=springCTF'*Mespring*springCTF;%整体坐标系下索的刚度矩阵
% 
% % 组集总体刚度矩阵和总体质量矩阵
% K(1:12,1:12)=K(1:12,1:12)+Kebeam;
% K(7:9,7:9)=K(7:9,7:9)+Kecable(1:3,1:3);
% K(7:9,13:15)=K(7:9,13:15)+Kecable(1:3,4:6);
% K(13:15,7:9)=K(13:15,7:9)+Kecable(4:6,1:3);
% K(13:15,13:15)=K(13:15,13:15)+Kecable(4:6,4:6);
% K(13:18,13:18)=K(13:18,13:18)+Kespring;%两个单元的矩阵
% 
% M(1:12,1:12)=M(1:12,1:12)+Mebeam;
% M(7:9,7:9)=M(7:9,7:9)+Mecable(1:3,1:3);
% M(7:9,13:15)=M(7:9,13:15)+Mecable(1:3,4:6);
% M(13:15,7:9)=M(13:15,7:9)+Mecable(4:6,1:3);
% M(13:15,13:15)=M(13:15,13:15)+Mecable(4:6,4:6);
% M(13:18,13:18)=M(13:18,13:18)+Mespring;%两个单元的矩阵
% 
% C = c1*M+c2*K;
% 
% Ka=K(7:15,7:15);
% Ma=M(7:15,7:15);
% Ca=C(7:15,7:15);
% 
% F = @(t) force(t,sumNode,everynodeDOF,sumBC_Displacement);
% [U,UD,UDD]=newMark(Ma,Ka,Ca,F,deltaT,T,delta,alpha);
% 
% % 绘图
% for i=1:9
%     dn = U(i,:);    
%     dV = UD(i,:);    
%     dA = UDD(i,:);  
%     tt = 0:deltaT:T;
%   %  hold on
%     figure(i)
%     subplot(2,2,1);
%     plot(tt,dn,'r','linewidth',0.1)
%     title('位移/转角')
%     subplot(2,2,2);
%     plot(tt,dV,'b','linewidth',0.1)
%     title('速度')
%     subplot(2,2,3);
%     plot(tt,dA,'g','linewidth',0.1)
%     title('加速度')
%    % hold off
% end
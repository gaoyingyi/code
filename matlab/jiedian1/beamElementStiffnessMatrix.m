%两节点十二自由度圆形梁单元，考虑剪切变形，考虑初应力，测试成功
%输入 内压力T是一个函数
%输入 L：是原长
%输出12×12刚度阵、质量阵
function [k,m]=beamElementStiffnessMatrix(T,L)
%% 参数
E=350000;                                                                   %弹性模量mpa
mu=0.3;                                                                     %泊松比
rho=1.5e-9;                                                                 %密度
r=5;                                                                        %圆截面半径为5mm
G=E/(2*(1+mu));                                                             %计算剪切模量
D=2*r;                                                                      %计算截面直径
A=pi*r^2;                                                                   %计算截面面积
Iy=pi*D^4/64;                                                               %计算圆截面惯性矩
Iz=pi*D^4/64;                                                               %计算圆截面惯性矩
Ip=pi*D^4/32;                                                               %计算圆截面极惯性矩
Ay=A/10*9;                                                                  %y方向有效抗剪面积
Az=A/10*9;                                                                  %Z方向有效抗剪面积
phaY=12*E*Iz/(G*Ay*L^2);                                                    %沿单元坐标y方向的横向力剪切影响系数
phaZ=12*E*Iy/(G*Az*L^2);                                                    %沿单元坐标z方向的横向力剪切影响系数
ry=sqrt(Iy/A);
rz=sqrt(Iz/A);
%% 刚度阵
%  线性与非线性有限元及其应用
%弹性刚度矩阵
kl=zeros(12);
k1=A*E/L;
k2=12*E*Iz/((1+phaY)*L^3);
k3=12*E*Iy/((1+phaZ)*L^3);
k4=G*Ip/L;
k5=(4+phaZ)*E*Iy/((1+phaZ)*L);
k6=(4+phaY)*E*Iz/((1+phaY)*L);
k7=-6*E*Iy/((1+phaZ)*L^2);
k8=6*E*Iz/((1+phaY)*L^2);
k9=(2-phaZ)*E*Iy/((1+phaZ)*L);
k10=(2-phaY)*E*Iz/((1+phaY)*L);
kl(1,1)=k1;kl(2,2)=k2;kl(3,3)=k3;kl(4,4)=k4;kl(5,5)=k5;kl(6,6)=k6;
kl(7,7)=k1;kl(8,8)=k2;kl(9,9)=k3;kl(10,10)=k4;kl(11,11)=k5;kl(12,12)=k6;
kl(5,3)=k7;kl(6,2)=k8;kl(7,1)=-k1;kl(8,2)=-k2;kl(8,6)=-k8;kl(9,3)=-k3;
kl(9,5)=-k7;kl(10,4)=-k4;kl(11,3)=k7;kl(11,5)=k9;kl(11,9)=-k7;kl(12,2)=k8;
kl(12,6)=k10;kl(12,8)=-k8;
kl(3,5)=k7;kl(2,6)=k8;kl(1,7)=-k1;kl(2,8)=-k2;kl(6,8)=-k8;kl(3,9)=-k3;
kl(5,9)=-k7;kl(4,10)=-k4;kl(3,11)=k7;kl(5,11)=k9;kl(9,11)=-k7;
kl(2,12)=k8;kl(6,12)=k10;kl(8,12)=-k8;                                     
%几何刚度矩阵，由应力状态引起的
g1=[0,0,      0;
    0,6/(5*L),0;
    0,0,      6/(5*L);];
g2=[0, 0,   0;
    0, 0,   L/10;
    0,-L/10,0;];
g3=[0,2*L/15,0;
    0,0,     2*L/15;
    0,0,     0;];
g4=[0,0,     0;
    0,-L/30, 0;
    0,0,    -L/30;];
kSigma=T.*[g1,   g2, -g1,  g2;
           g2',  g3,  g2,  g4;
          -g1',  g2', g1, -g2;
           g2',  g4',-g2', g3;];                                          
k=kl+kSigma;
%% 一致质量阵
m=zeros(12);
me=rho*A*L;
Ay=(13/35+7/10*phaZ+1/3*phaZ^2+6/5*(ry/L)^2)/(1+phaZ)^2;
Az=(13/35+7/10*phaY+1/3*phaY^2+6/5*(rz/L)^2)/(1+phaY)^2;
By=(9/70+3/10*phaZ+1/6*phaZ^2-6/5*(ry/L)^2)/(1+phaZ)^2;
Bz=(9/70+3/10*phaY+1/6*phaY^2-6/5*(rz/L)^2)/(1+phaY)^2;
Cy=(11/210+11/120*phaZ+1/24*phaZ^2+(1/10-1/2*phaZ)*(ry/L)^2)*L/(1+phaZ)^2;
Cz=(11/210+11/120*phaY+1/24*phaY^2+(1/10-1/2*phaY)*(rz/L)^2)*L/(1+phaY)^2;
Dy=(13/420+3/40*phaZ+1/24*phaZ^2-(1/10-1/2*phaZ)*(ry/L)^2)*L/(1+phaZ)^2;
Dz=(13/420+3/40*phaY+1/24*phaY^2-(1/10-1/2*phaY)*(rz/L)^2)*L/(1+phaY)^2;
Ey=(1/105+1/60*phaZ+1/120*phaZ^2+(2/15+1/6*phaZ+1/3*phaZ^2)*(ry/L)^2)*L^2/(1+phaZ)^2;
Ez=(1/105+1/60*phaY+1/120*phaY^2+(2/15+1/6*phaY+1/3*phaY^2)*(rz/L)^2)*L^2/(1+phaY)^2;
Fy=-(1/140+1/60*phaZ+1/120*phaZ^2+(1/30+1/6*phaZ-1/6*phaZ^2)*(ry/L)^2)*L^2/(1+phaZ)^2;
Fz=-(1/140+1/60*phaY+1/120*phaY^2+(1/30+1/6*phaY-1/6*phaY^2)*(rz/L)^2)*L^2/(1+phaY)^2; 
m(1,1)=1/3*me;m(2,2)=Az*me;m(3,3)=Ay*me;m(4,4)=1/3*Ip/A*me;m(5,5)=Ey*me;m(6,6)=Ez*me;m(7,7)=1/3*me;m(8,8)=Az*me;m(9,9)=Ay*me;m(10,10)=1/3*Ip/A*me;m(11,11)=Ey*me;m(12,12)=Ez*me;
m(5,3)=-Cy*me;m(6,2)=Cz*me;m(11,9)=Cy*me;m(12,8)=-Cz*me;m(7,1)=1/6*me;m(8,2)=Bz*me;m(9,3)=By*me;m(10,4)=1/6*Ip/A*me;m(11,5)=Fy*me;m(12,6)=Fz*me;m(12,2)=-Dz*me;m(8,6)=-m(12,2);m(11,3)=Dy*me;m(9,5)=-m(11,3);
m(3,5)=-Cy*me;m(2,6)=Cz*me;m(9,11)=Cy*me;m(8,12)=-Cz*me;m(1,7)=1/6*me;m(2,8)=Bz*me;m(3,9)=By*me;m(4,10)=1/6*Ip/A*me;m(5,11)=Fy*me;m(6,12)=Fz*me;m(2,12)=-Dz*me;m(6,8)=-m(12,2);m(3,11)=Dy*me;m(5,9)=-m(11,3);
end
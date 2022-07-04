clc;clear;close all;
Le=4;
d=10;
v=0.3;
E=350000;
p=1.5e-9;
Iy=pi/64*d^4;
Iz=pi/64*d^4;
Jk=Iy+Iz;                                        %梁极惯性距
A=0.25*d*d*pi;
G=E/2/(1+v);  
ry=sqrt(Iy/A);
rz=sqrt(Iz/A);
pha_y=12*E*Iz/G/(0.9*A)/Le/Le;                      %沿单元坐标y方向的横向力剪切影响系数
pha_z=12*E*Iy/G/(0.9*A)/Le/Le;                      %沿单元坐标z方向的横向力剪切影响系数
Ay=(13/35+7/10*pha_z+1/3*pha_z^2+6/5*(ry/Le)^2)/(1+pha_z)^2;
Az=(13/35+7/10*pha_y+1/3*pha_y^2+6/5*(rz/Le)^2)/(1+pha_y)^2;
By=(9/70+3/10*pha_z+1/6*pha_z^2-6/5*(ry/Le)^2)/(1+pha_z)^2;
Bz=(9/70+3/10*pha_y+1/6*pha_y^2-6/5*(rz/Le)^2)/(1+pha_y)^2;
Cy=(11/210+11/120*pha_z+1/24*pha_z^2+(1/10-1/2*pha_z)*(ry/Le)^2)*Le/(1+pha_z)^2;
Cz=(11/210+11/120*pha_y+1/24*pha_y^2+(1/10-1/2*pha_y)*(rz/Le)^2)*Le/(1+pha_y)^2;
Dy=(13/420+3/40*pha_z+1/24*pha_z^2-(1/10-1/2*pha_z)*(ry/Le)^2)*Le/(1+pha_z)^2;
Dz=(13/420+3/40*pha_y+1/24*pha_y^2-(1/10-1/2*pha_y)*(rz/Le)^2)*Le/(1+pha_y)^2;
Ey=(1/105+1/60*pha_z+1/120*pha_z^2-(2/15+1/6*pha_z+1/3*pha_z^2)*(ry/Le)^2)*Le^2/(1+pha_z)^2;
Ez=(1/105+1/60*pha_y+1/120*pha_y^2-(2/15+1/6*pha_y+1/3*pha_y^2)*(rz/Le)^2)*Le^2/(1+pha_y)^2;
Fy=(1/140+1/60*pha_z+1/120*pha_z^2+(1/30+1/6*pha_z-1/6*pha_z^2)*(ry/Le)^2)*Le^2/(1+pha_z)^2;
Fz=(1/140+1/60*pha_y+1/120*pha_y^2+(1/30+1/6*pha_y-1/6*pha_y^2)*(rz/Le)^2)*Le^2/(1+pha_y)^2;
me=p*A*Le;  
%单元质量矩阵中的共同系数   
Mii(1,1)=1/3*me;
Mii(2,2)=Az*me;
Mii(3,3)=Ay*me;
Mii(4,4)=1/3*Jk/A*me;
Mii(5,5)=2*Ey*me;
Mii(6,6)=2*Ez*me;
Mii(6,2)=-1/2*Cz*me;
Mii(2,6)=Mii(6,2);
Mii(5,3)=1/2*Cy*me;
Mii(3,5)=Mii(5,3);

Mjj(1,1)=1/3*me;
Mjj(2,2)=Az*me;
Mjj(3,3)=Ay*me;
Mjj(4,4)=1/3*Jk/A*me;
Mjj(5,5)=2*Ey*me;
Mjj(6,6)=2*Ez*me;
Mjj(6,2)=1/2*Cz*me;
Mjj(2,6)=Mjj(6,2);
Mjj(5,3)=-1/2*Cy*me;
Mjj(3,5)=Mjj(5,3);

Mji(1,1)=1/6*me;
Mji(2,2)=Bz*me;
Mji(3,3)=By*me;
Mji(4,4)=-1/6*Jk/A*me/2;
Mji(5,5)=Fy*me;
Mji(6,6)=Fz*me;
Mji(6,2)=-1/2*Dz*me;
Mji(2,6)=-Mji(6,2);
Mji(5,3)=1/2*Dy*me;
Mji(3,5)=-Mji(5,3);

Mij=(Mji)';

M=[Mii,Mij;Mji,Mjj];
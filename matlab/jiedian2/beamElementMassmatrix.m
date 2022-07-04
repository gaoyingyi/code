%% 输出：梁单元一致质量阵
function m=beamElementMassmatrix(L)
%% 参数
global EBeam muBeam rhoBeam ;
G=EBeam/(2*(1+muBeam));                                                     %计算剪切模量
r=5;                                                                        %横截面半径
D=2*r;                                                                      %横截面直径
A=pi*r^2;                                                                   %横截面面积
Iy=pi*D^4/64;                                                               %横截面y惯性矩
Iz=pi*D^4/64;                                                               %横截面z惯性矩
Ip=pi*D^4/32;                                                               %横截面极惯性矩
Ay=A/10*9;                                                                  %y方向有效抗剪面积
Az=A/10*9;                                                                  %Z方向有效抗剪面积
phaY=12*EBeam*Iz/(G*Ay*L^2);                                                %y方向剪切影响系数
phaZ=12*EBeam*Iy/(G*Az*L^2);                                                %z方向剪切影响系数
ry=sqrt(Iy/A);
rz=sqrt(Iz/A);
%-------------------------------------------------------------------------------------------
me=rhoBeam*A*L;
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
%梁单元一致质量阵，经过检验
m=zeros(12);
m(1,1)=1/3*me;m(2,2)=Az*me;m(3,3)=Ay*me;m(4,4)=1/3*Ip/A*me;m(5,5)=Ey*me;m(6,6)=Ez*me;
m(7,7)=1/3*me;m(8,8)=Az*me;m(9,9)=Ay*me;m(10,10)=1/3*Ip/A*me;m(11,11)=Ey*me;m(12,12)=Ez*me;
m(5,3)=-Cy*me;m(6,2)=Cz*me;m(11,9)=Cy*me;m(12,8)=-Cz*me;m(7,1)=1/6*me;m(8,2)=Bz*me;m(9,3)=By*me;
m(10,4)=1/6*Ip/A*me;m(11,5)=Fy*me;m(12,6)=Fz*me;m(12,2)=-Dz*me;m(8,6)=-m(12,2);
m(11,3)=Dy*me;m(9,5)=-m(11,3);
m(3,5)=-Cy*me;m(2,6)=Cz*me;m(9,11)=Cy*me;m(8,12)=-Cz*me;m(1,7)=1/6*me;m(2,8)=Bz*me;m(3,9)=By*me;
m(4,10)=1/6*Ip/A*me;m(5,11)=Fy*me;m(6,12)=Fz*me;m(2,12)=-Dz*me;m(6,8)=-m(12,2);
m(3,11)=Dy*me;m(5,9)=-m(11,3);
end
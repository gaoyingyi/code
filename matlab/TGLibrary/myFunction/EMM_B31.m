%输出：梁单元一致质量阵
%杜柏松,项海帆,葛耀君,朱乐东.剪切效应梁单元刚度和质量矩阵的推导及应用[J].重庆交通大学学报(自然科学版),2008(04):502-507.
function m=EMM_B31(E,mu,rho,r,nodeCoordinateMatrix)

xi=nodeCoordinateMatrix(1,1);
yi=nodeCoordinateMatrix(1,2);
zi=nodeCoordinateMatrix(1,3);
xj=nodeCoordinateMatrix(2,1);
yj=nodeCoordinateMatrix(2,2);
zj=nodeCoordinateMatrix(2,3);
L=sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);%梁单元的长度
%% 参数
G=E/(2*(1+mu));%计算剪切模量
D=2*r;%横截面直径
A=pi*r^2;%横截面面积
Iy=pi*D^4/64;%横截面y惯性矩
Iz=pi*D^4/64;%横截面z惯性矩
J=pi*D^4/32;%横截面极惯性矩
Ay=A/10*9;%y方向有效抗剪面积
Az=A/10*9;%Z方向有效抗剪面积
phaY=12*E*Iz/(G*Ay*L^2);%y方向剪切影响系数
phaZ=12*E*Iy/(G*Az*L^2);%z方向剪切影响系数
ry=sqrt(Iy/A);
rz=sqrt(Iz/A);
%-------------------------------------------------------------------------------------------
me=rho*A*L;
Ay=(13/35+7/10*phaZ+1/3*phaZ^2+6/5*(ry/L)^2)/(1+phaZ)^2;
Az=(13/35+7/10*phaY+1/3*phaY^2+6/5*(rz/L)^2)/(1+phaY)^2;
By=(9/70+3/10*phaZ+1/6*phaZ^2-6/5*(ry/L)^2)/(1+phaZ)^2;
Bz=(9/70+3/10*phaY+1/6*phaY^2-6/5*(rz/L)^2)/(1+phaY)^2;
Cy=(11/210+11/120*phaZ+1/24*phaZ^2+(1/10-1/2*phaZ)*(ry/L)^2)*L/(1+phaZ)^2;
Cz=(11/210+11/120*phaY+1/24*phaY^2+(1/10-1/2*phaY)*(rz/L)^2)*L/(1+phaY)^2;
Dy=(13/420+3/40*phaZ+1/24*phaZ^2-(1/10-1/2*phaZ)*(ry/L)^2)*L/(1+phaZ)^2;
Dz=(13/420+3/40*phaY+1/24*phaY^2-(1/10-1/2*phaY)*(rz/L)^2)*L/(1+phaY)^2;
Ey=(1/105+1/60*phaZ+1/120*phaZ^2+(2/15+1/6*phaZ+1/3*phaZ^2)*(ry/L)^2)*L^2/...
(1+phaZ)^2;
Ez=(1/105+1/60*phaY+1/120*phaY^2+(2/15+1/6*phaY+1/3*phaY^2)*(rz/L)^2)*L^2/...
(1+phaY)^2;
Fy=-(1/140+1/60*phaZ+1/120*phaZ^2+(1/30+1/6*phaZ-1/6*phaZ^2)*(ry/L)^2)*...
L^2/(1+phaZ)^2;
Fz=-(1/140+1/60*phaY+1/120*phaY^2+(1/30+1/6*phaY-1/6*phaY^2)*(rz/L)^2)*...
L^2/(1+phaY)^2;
%梁单元一致质量阵，经过检验
m=zeros(12);
m(1,1)=1/3;m(2,2)=Az;m(3,3)=Ay;  m(4,4)=1/3*J/A;  m(5,5)=Ey;  m(6,6)=Ez;
m(7,7)=1/3;m(8,8)=Az;m(9,9)=Ay;m(10,10)=1/3*J/A;m(11,11)=Ey;m(12,12)=Ez;
m(5,3)=-Cy;m(6,2)=Cz;m(7,1)=1/6;m(8,2)=Bz;m(8,6)=-m(12,2);m(9,3)=By;m(9,5)=-m(11,3);m(10,4)=1/6*J/A;m(11,3)=Dy;m(11,5)=Fy;m(11,9)=Cy;m(12,2)=-Dz;m(12,6)=Fz;m(12,8)=-Cz;
m(3,5)=-Cy;m(2,6)=Cz;m(1,7)=1/6;m(2,8)=Bz;m(6,8)=-m(12,2);m(3,9)=By;m(5,9)=-m(11,3);m(4,10)=1/6*J/A;m(3,11)=Dy;m(5,11)=Fy;m(9,11)=Cy;m(2,12)=-Dz;m(6,12)=Fz;m(8,12)=-Cz;
m=me*m;

%% ========================================================================
%局部到全局坐标变换                                                                    
ll=(xj-xi)/L;
mm=(yj-yi)/L;
nn=(zj-zi)/L;
DD=sqrt(ll^2+mm^2);
if nn==1                           %判断是否和全局Z轴方向相同
    lambda=[0,0,1;
            0,1,0;
            -1,0,0;];
    elseif nn==-1                   %判断是否和全局Z轴方向相反
    lambda=[0,0,-1;
            0,1,0;
            1,0,0;];
    else                            %一般情况
    lambda=[     ll,       mm,  mm;
              -mm/DD,    -ll/DD,   0;
           -ll*nn/DD, -mm*nn/DD,   DD;];
end
T=[lambda,zeros(3),zeros(3),zeros(3);
   zeros(3),lambda,zeros(3),zeros(3);
   zeros(3),zeros(3),lambda,zeros(3);
   zeros(3),zeros(3),zeros(3),lambda;]; 
m=T'*m*T;
end
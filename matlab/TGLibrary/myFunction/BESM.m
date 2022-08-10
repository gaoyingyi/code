%功能：给定梁的材料参数和在全局坐标系下的坐标，局部坐标系下
function k=BESM(E,mu,nodeCoordinateMatrix,nodeDisplacementMatrix)
% E:弹性模量
% mu 泊松比
% nodeCoordinateMatrix全局坐标系下节点坐标矩阵2×3，存储x_i，y_i，z_i，x_j，y_j，z_j
% nodeDisplacementMatrix局部节点位移矩阵12×1，存储u_i，v_i，w_i，θ_xi，θ_yi，θ_zi，
%                                        u_j，v_j，w_j，θ_xj，θ_yj，θ_zj
%% ========================================================================
u_i=nodeDisplacementMatrix(1,1);
v_i=nodeDisplacementMatrix(2,1);
w_i=nodeDisplacementMatrix(3,1);
theta_xi=nodeDisplacementMatrix(4,1);
theta_yi=nodeDisplacementMatrix(5,1);
theta_zi=nodeDisplacementMatrix(6,1);
u_j=nodeDisplacementMatrix(7,1);
v_j=nodeDisplacementMatrix(8,1);
w_j=nodeDisplacementMatrix(9,1);
theta_xj=nodeDisplacementMatrix(10,1);
theta_yj=nodeDisplacementMatrix(11,1);
theta_zj=nodeDisplacementMatrix(12,1);

x_i=nodeCoordinateMatrix(1,1);
y_i=nodeCoordinateMatrix(1,2);
z_i=nodeCoordinateMatrix(1,3);
x_j=nodeCoordinateMatrix(2,1);
y_j=nodeCoordinateMatrix(2,2);
z_j=nodeCoordinateMatrix(2,3);
l=sqrt((x_i-x_j)^2+(y_i-y_j)^2+(z_i-z_j)^2);%梁的长度
%% ========================================================================
%线性刚度阵预备参数
G=E/(2*(1+mu));%计算剪切模量
%--------------------------截面属性-----------------------------------------
r=5;%圆截面半径
D=2*r;%计算截面直径
A=pi*r^2;%计算截面面积
Iy=pi*D^4/64;%计算圆截面惯性矩
Iz=pi*D^4/64;%计算圆截面惯性矩
J=pi*D^4/32;%计算圆截面极惯性矩
Ay=A/10*9;%y方向有效抗剪面积
Az=A/10*9;%Z方向有效抗剪面积

phaY=12*E*Iz/(G*Ay*l^2);%沿单元坐标y方向的横向力剪切影响系数
phaZ=12*E*Iy/(G*Az*l^2);%沿单元坐标z方向的横向力剪切影响系数
H=2;
%% ========================================================================
%线性刚度矩阵
k0=zeros(12);%初始化线性刚度阵
k0(1,1)=A*E/l;
k0(2,2)=12*E*Iz/((1+phaY)*l^3);
k0(3,3)=12*E*Iy/((1+phaZ)*l^3);
k0(4,4)=G*J/l;
k0(5,5)=(4+phaZ)*E*Iy/((1+phaZ)*l);
k0(6,6)=(4+phaY)*E*Iz/((1+phaY)*l);
k0(5,3)=-6*E*Iy/((1+phaZ)*l^2);
k0(6,2)=6*E*Iz/((1+phaY)*l^2);
k0(11,5)=(2-phaZ)*E*Iy/((1+phaZ)*l);
k0(12,6)=(2-phaY)*E*Iz/((1+phaY)*l);%这里列出了矩阵中10个不同的值，其他的值用这个赋予。
%下三角矩阵的其他值
k0(7,7)=k0(1,1);k0(8,8)=k0(2,2);k0(9,9)=k0(3,3);k0(10,10)=k0(4,4);k0(11,11)=k0(5,5);k0(12,12)=k0(6,6);
k0(7,1)=-k0(1,1);k0(8,2)=-k0(2,2);k0(8,6)=-k0(6,2);k0(9,3)=-k0(3,3);k0(9,5)=-k0(5,3);k0(10,4)=-k0(4,4);k0(11,3)=k0(5,3);k0(11,9)=-k0(5,3);k0(12,2)=k0(6,2);k0(12,8)=-k0(6,2);
% 对称矩阵，给矩阵上三角赋值
k0(1,7)=-k0(1,1);k0(2,8)=-k0(2,2);k0(6,8)=-k0(6,2);k0(3,9)=-k0(3,3);k0(5,9)=-k0(5,3);k0(4,10)=-k0(4,4);k0(3,11)=k0(5,3);k0(9,11)=-k0(5,3);k0(2,12)=k0(6,2);k0(8,12)=-k0(6,2);
k0(3,5)=k0(5,3);k0(2,6)=k0(6,2);k0(5,11)=k0(11,5);k0(6,12)=k0(12,6);
%% ========================================================================
%大位移刚度矩阵
kL=zeros(12);%初始化大位移刚度阵
kL(2,1)=6*E*A*(v_i-v_j)/(5*l^2)-E*A*(theta_zi+theta_zj)/(10*l);%
kL(3,1)=6*E*A*(w_i-w_j)/(5*l^2)+E*A*(theta_yi+theta_yj)/(10*l);% 
kL(4,1)=E*J*(theta_xi-theta_xj)/l^2;%
kL(5,1)=E*A*(w_i-w_j)/(10*l)-2*E*A*theta_yi/15+E*A*theta_yj/30;%
kL(6,1)=E*A*(v_i-v_j)/(10*l)-2*E*A*theta_zi/15+E*A*theta_zj/30;%
kL(11,1)=E*A*(w_i-w_j)/(10*l)+E*A*theta_yi/30-2*E*A*theta_yj/15;%
kL(12,1)=-E*A*(v_i-v_j)/(10*l)+E*A*theta_zi/30-2*E*A*theta_zj/15;%
kL(2,2)=72*E*A/(35*l^3)*(v_i-v_j)^2+18*E*A/(35*l^2)*(v_i-v_j)*(theta_zi+theta_zj)+3*E*A/(35*l)*(theta_zi^2+theta_zj^2)+(6*E*Iz)/(5*l)*(theta_xi-theta_xj)^2;%
kL(3,2)=72*E*A/(35*l^3)*(v_i-v_j)*(w_i-w_j)-9*E*A/(35*l^2)*(v_i-v_j)*(theta_yi+theta_yj)+9*E*A/(35*l^2)*(w_i-w_j)*(theta_zi+theta_zj)-3*E*A/(35*l)*(theta_zi*theta_yi+theta_zj*theta_yj);%
kL(4,2)=E*(Iy+2*Iz)/(10*l^3)*(theta_xi-theta_xj)*(12*(v_i-v_j)+l*(theta_zi+theta_zj));%
kL(5,2)=-3*E*A/(35*l)*theta_zi*(w_i-w_j)+3*E*A/(35*l)*theta_yi*(v_i-v_j)-E*A/140*theta_zi*(theta_yi-theta_yj)+E*A/140*theta_zj*(theta_yi+theta_yj)-9*E*A/(35*l^2)*(v_i-v_j)*(w_i-w_j);%
kL(6,2)=9*E*A/(35*l^2)*(v_i-v_j)^2+6*E*A/(35*l)*theta_zi*(v_i-v_j)-E*A/140*(theta_zi^2-theta_zj^2)+E*A/70*theta_zi*theta_zj+(E*Iz)/(10*l^2)*(theta_xi-theta_xj)^2;%
kL(11,2)=-9*E*A/(35*l^2)*(v_i-v_j)*(w_i-w_j)-3*E*A/(35*l)*theta_zj*(w_i-w_j)+3*E*A/(35*l)*theta_yj*(v_i-v_j)+E*A/140*theta_xj*(theta_yi-theta_yj)+E*A/140*theta_zi*(theta_yi+theta_yj);%
kL(12,2)=9*E*A/(35*l^2)*(v_i-v_j)^2+6*E*A/(35*l)*theta_zj*(v_i-v_j)+E*A/140*(theta_zi^2-theta_xj^2)+E*A/70*theta_zi*theta_zj+(E*Iz)/(10*l^2)*(theta_xi-theta_xj)^2;%
kL(3,3)=72*E*A/(35*l^3)*(w_i-w_j)^2-18*E*A/(35*l^2)*(w_i-w_j)*(theta_yi+theta_yj)+3*E*A/(35*l)*(theta_yi^2+theta_yj^2)+(42*E*Ix)/(35*l^3)*(theta_xi-theta_xj)^2;%
kL(4,3)=E(2*Iy+Iz)/(10*l^3)*(theta_xi-theta_xj)*(12*(w_i-w_j)-l*(theta_yi+theta_yj));%
kL(5,3)=-9*E*A/(35*l^2)*(w_i-w_j)^2+6*E*A/(35*l)*theta_yi*(w_i-w_j)+E*A/140*(theta_yi^2-theta_yj^2)-E*A/70*theta_yi*theta_yj-(E*Iy)/(10*l^2)*(theta_xi-theta_xj)^2;%
kL(6,3)=9*E*A/(35*l^2)*(v_i-v_j)*(w_i-w_j)+3*E*A/(35*l)*theta_zi*(w_i-w_j)-3*E*A/(35*l)*theta_yi*(v_i-v_j)+E*A/140*theta_zi*(theta_yi-theta_yj)-E*A/140*theta_zj*(theta_yi+theta_yj);%
kL(11,3)=-9*E*A/35*l^2*(w_i-w_j)^2+6*E*A/(35*l)*theta_yj*(w_i-w_j)-E*A/140*(theta_yi^2-theta_yj^2)-E*A/70*theta_yi*theta_yj-(E*Iy)/(10*l^2)*(theta_xi-theta_xj)^2;%
kL(12,3)=9*E*A/(35*l^2)*(v_i-v_j)*(w_i-w_j)+3*E*A/(35*l)*theta_zj*(w_i-w_j)-3*E*A/(35*l)*theta_yi (v_i-v_j)-E*A/140*theta_zj*(theta_yi-theta_yj)-E*A/140*theta_zi*(theta_yi+theta_yj);%
kL(4,4)=(6*E*Iz)/(5*l^3)*(v_i-v_j)^2+(6*E*Iy)/(5*l^3)*(w_i-w_j)^2+(E*A^3)/(72*l^3)*(theta_xi-theta_xj)^2+(E*A*W^4)/(80*l^3)*(theta_xi-theta_xj)^2+(E*A*H^4)/(80*l^3)*(theta_xi-theta_xj)^2+(E*Iz)/(5*l^2)*(v_i-v_j)*(theta_zi+theta_zj)-(E*Iy)/(5*l^2)*(w_i-w_j)*(theta_yi+theta_yj)-(E*Iz)/15*theta_zi*theta_zj-(E*Iy)/(15*l)*theta_yi*theta_yj+(2*E*Iy)/(15*l)*(theta_yi^2+theta_yj^2)+(2*E*Iz)/(15*l)*(theta_zi^2+theta_zj^2);%
kL(5,4)=E*(2*Iy+Iz)/(30*l^2)*(theta_xi-theta_xj)*(-3*(w_i-w_j)+l*(4*theta_yi-theta_yj));%
kL(6,4)=E*(Iy+2*Iz)/(30*l^2)*(theta_xi-theta_xj)*( 3*(v_i-v_j)+l*(4*theta_zi-theta_zj));%
kL(11,4)=-E*(2*Iy+Iz)/(30*l^2)*(theta_xi-theta_xj)*(3*(w_i-w_j)+l*(theta_yi-4*theta_yj));%
kL(12,4)= E*(Iy+2*Iz)/(30*l^2)*(theta_xi-theta_xj)*(3*(v_i-v_j)-l*(theta_zi-4*theta_zj));%
kL(5,5)=  3*E*A/(35*l)*(w_i-w_j)^2+E*A/70*(w_i-w_j)*(theta_yi-theta_yj)+2*E*A*l/35*theta_yi^2+E*A*l/210*theta_yj^2-E*A*l/70*theta_yi*theta_yj+(2*E*Iy)/(15*l)*(theta_xi-theta_xj)^2;%
kL(6,5)=-3*E*A/35*l (v_i-v_j)*(w_i-w_j)-E*A/140*(v_i-v_j)*(theta_yi-theta_yj)+E*A/140*(w_i-w_j)*(theta_zi-theta_zj)-E*A*l/140*(theta_zi*theta_yj+theta_zj*theta_yi)+2*E*A*l/35*theta_zi*theta_yi+E*A*l/210*theta_zj*theta_yj;%
kL(11,5)=-E*A/70*(w_i-w_j)*(theta_yi+theta_yj)-E*A*l/140*(theta_yi^2+theta_yj^2)+E*A*l/105*theta_yi*theta_yj-(E*Iy)/(30*l)*(theta_xi-theta_xj)^2;%
kL(12,5)=-E*A/140*(w_i-w_j)*(theta_zi+theta_zj)+E*A/140*(v_i-v_j)*(theta_yi+theta_yj)-E*A*l/140*(theta_zi*theta_yi+theta_zj*theta_yj)+E*A*l/210*(theta_zi*theta_yj+theta_zj*theta_yj);%
kL(6,6)=3*E*A/(35*l)*(v_i-v_j)^2-E*A/70*(v_i-v_j)*(theta_zi-theta_zj)+2*E*A*l/35*theta_zi^2+E*A*l/210*theta_zj^2-E*A*l/70*theta_zi*theta_zj+(2*E*Iz)/(15*l)*(theta_xi-theta_xj)^2;%
kL(11,6)=-E*A/140*(w_i-w_j)*(theta_zi+theta_zj)+E*A/140*(v_i-v_j)*(theta_yi+theta_yj)-E*A*l/140*(theta_zi*theta_yi+theta_zj*theta_yj)+E*A*l/210*(theta_zi*theta_yj+theta_zj*theta_yi);%
kL(12,6)= E*A/70*(v_i-v_j)*(theta_zi+theta_zj)-E*A*l/140*(theta_zi^2+theta_zj^2)+E*A*l/105*theta_zi*theta_zj-(E*Iz)/(30*l)*(theta_xi-theta_xj)^2;%
kL(11,11)=3*E*A/(35*l)*(w_i-w_j)^2-E*A/70*(w_i-w_j)*(theta_yi-theta_yj)+2*E*A*l/35*theta_yj^2+E*A*l/210*theta_yi^2-E*A*l/70*theta_yi*theta_yj+(2*E*Iy)/(15*l)*(theta_xi-theta_xj)^2;%
kL(12,11)=-E*A/140*(w_i-w_j)*(theta_zi-theta_zj)+E*A/140*(v_i-v_j)*(theta_yi-theta_yj)-3*E*A/(35*l)*(w_i-w_j)*(v_i-v_j)-E*A*l/140*(theta_zi*theta_yj+theta_zj*theta_yi)+E*A*l/210*theta_zi*theta_yi+2*E*A*l/35*theta_zj*theta_yj;%
kL(12,12)=3*E*A/(35*l)*(v_i-v_j)^2+E*A/70*(v_i-v_j)*(theta_zi-theta_zj)+E*A*l/210*theta_zi^2+2*E*A*l/35*theta_zj^2-E*A*l/70*theta_zi*theta_zj+(2*E*Iz)/(15*l)*(theta_xi-theta_xj)^2;%
kL(8,8)=kL(2,2);kL(9,9)=kL(3,3);kL(10,10)=kL(4,4);%主对角线赋值
kL(1,2)=kL(2,1);kL(1,3)=kL(3,1);kL(1,4)=kL(4,1);kL(1,5)=kL(5,1);kL(1,6)=kL(6,1);%
kL(2,3)=kL(3,2);kL(2,4)=kL(4,2);kL(2,5)=kL(5,2);kL(2,6)=kL(6,2);%
kL(3,4)=kL(4,3);kL(3,5)=kL(5,3);kL(3,6)=kL(6,3);%
kL(4,5)=kL(5,4);kL(4,6)=kL(6,4);%
kL(5,6)=kL(6,5);%块1赋值
kL(7,2)=-kL(2,1);kL(7,3)=-kL(3,1);kL(7,4)=-kL(4,1);kL(7,5)=-kL(5,1);kL(7,6)=-kL(6,1);%
kL(8,1)=-kL(2,1);kL(8,2)=-kL(2,2);kL(8,3)=-kL(3,2);kL(8,4)=-kL(4,2);kL(8,5)=-kL(5,2);kL(8,6)=-kL(6,2);%
kL(9,1)=-kL(3,1);kL(9,2)=-kL(3,2);kL(9,3)=-kL(3,3);kL(9,4)=-kL(4,3);kL(9,5)=-kL(5,3);kL(9,6)=-kL(6,3);%
kL(10,1)=-kL(4,1);kL(10,2)=-kL(4,2);kL(10,3)=-kL(4,3);kL(10,4)=-kL(4,4);kL(10,5)=-kL(4,5);kL(10,6)=-kL(4,6);%块2赋值
kL(8,7)=kL(2,1);kL(9,7)=kL(3,1);kL(9,8)=kL(3,2);kL(10,7)=kL(4,1);kL(10,8)=kL(4,2);kL(10,9)=kL(4,3);%块3赋值
kL(11,7)=-kL(11,1);kL(11,8)=-kL(11,2);kL(11,9)=-kL(11,3);kL(11,10)=-kL(11,4);%
kL(12,7)=-kL(12,1);kL(12,8)=-kL(12,2);kL(12,9)=-kL(12,3);kL(12,10)=-kL(12,4);%块4赋值
kL(1,2)=kL(2,1);kL(1,3)=kL(3,1);kL(1,4)=kL(4,1);kL(1,5)=kL(5,1);kL(1,6)=kL(6,1);kL(1,7)=kL(7,1);
kL(1,8)=kL(8,1);kL(1,9)=kL(9,1);kL(1,10)=kL(10,1);kL(1,11)=kL(11,1);kL(1,12)=kL(12,1);kL(2,3)=kL(3,2);
kL(2,4)=kL(4,2);kL(2,5)=kL(5,2);kL(2,6)=kL(6,2);kL(2,7)=kL(7,2);kL(2,8)=kL(8,2);kL(2,9)=kL(9,2);kL(2,10)=kL(10,2);
kL(2,11)=kL(11,2);kL(2,12)=kL(12,2);kL(3,4)=kL(4,3);kL(3,5)=kL(5,3);kL(3,6)=kL(6,3);kL(3,7)=kL(7,3);
kL(3,8)=kL(8,3);kL(3,9)=kL(9,3);kL(3,10)=kL(10,3);kL(3,11)=kL(11,3);kL(3,12)=kL(12,3);kL(4,5)=kL(5,4);
kL(4,6)=kL(6,4);kL(4,7)=kL(7,4);kL(4,8)=kL(8,4);kL(4,9)=kL(9,4);kL(4,10)=kL(10,4);kL(4,11)=kL(11,4);
kL(4,12)=kL(12,4);kL(5,6)=kL(6,5);kL(5,7)=kL(7,5);kL(5,8)=kL(8,5);kL(5,9)=kL(9,5);kL(5,10)=kL(10,5);
kL(5,11)=kL(11,5);kL(5,12)=kL(12,5);kL(6,7)=kL(7,6);kL(6,8)=kL(8,6);kL(6,9)=kL(9,6);kL(6,10)=kL(10,6);
kL(6,11)=kL(11,6);kL(6,12)=kL(12,6);kL(7,8)=kL(8,7);kL(7,9)=kL(9,7);kL(7,10)=kL(10,7);kL(7,11)=kL(11,7);
kL(7,12)=kL(12,7);kL(8,9)=kL(9,8);kL(8,10)=kL(10,8);kL(8,11)=kL(11,8);kL(8,12)=kL(12,8);kL(9,10)=kL(10,9);
kL(9,11)=kL(11,9);kL(9,12)=kL(12,9);kL(10,11)=kL(11,10);kL(10,12)=kL(12,10);kL(11,12)=kL(12,11);%12×12上三角对称矩阵赋值
%% ========================================================================
%初应力刚度矩阵
kSigma0=zeros(12);%初始化初应力刚度阵
kSigma0(2,2)=6/(5*l);
kSigma0(4,4)=(Iy+Iz)/(l*A);
kSigma0(5,5)=2*l/15;
kSigma0(6,2)=1/10;
kSigma0(11,5)=1/30;
kSigma0(3,3)=kSigma0(2,2);kSigma0(6,6)=kSigma0(5,5);kSigma0(8,8)=kSigma0(2,2);kSigma0(9,9)=kSigma0(2,2);kSigma0(10,10)=kSigma0(4,4);kSigma0(11,11)=kSigma0(5,5);kSigma0(12,12)=kSigma0(5,5);%对角线赋值
kSigma0(5,3)=-kSigma0(6,2);kSigma0(8,2)=-kSigma0(2,2);kSigma0(8,6)=-kSigma0(6,2);kSigma0(9,3)=-kSigma0(6,2);kSigma0(9,5)=-kSigma0(2,2);
kSigma0(10,4)=-kSigma0(4,4);kSigma0(11,3)=-kSigma0(6,2);kSigma0(11,9)=kSigma0(6,2);kSigma0(12,2)=kSigma0(6,2);kSigma0(12,6)=-kSigma0(11,5);kSigma0(12,8)=-kSigma0(6,2);%下三角矩阵赋值
kSigma0(1,2)=kSigma0(2,1);kSigma0(1,3)=kSigma0(3,1);kSigma0(1,4)=kSigma0(4,1);kSigma0(1,5)=kSigma0(5,1);kSigma0(1,6)=kSigma0(6,1);kSigma0(1,7)=kSigma0(7,1);
kSigma0(1,8)=kSigma0(8,1);kSigma0(1,9)=kSigma0(9,1);kSigma0(1,10)=kSigma0(10,1);kSigma0(1,11)=kSigma0(11,1);kSigma0(1,12)=kSigma0(12,1);kSigma0(2,3)=kSigma0(3,2);
kSigma0(2,4)=kSigma0(4,2);kSigma0(2,5)=kSigma0(5,2);kSigma0(2,6)=kSigma0(6,2);kSigma0(2,7)=kSigma0(7,2);kSigma0(2,8)=kSigma0(8,2);kSigma0(2,9)=kSigma0(9,2);kSigma0(2,10)=kSigma0(10,2);
kSigma0(2,11)=kSigma0(11,2);kSigma0(2,12)=kSigma0(12,2);kSigma0(3,4)=kSigma0(4,3);kSigma0(3,5)=kSigma0(5,3);kSigma0(3,6)=kSigma0(6,3);kSigma0(3,7)=kSigma0(7,3);
kSigma0(3,8)=kSigma0(8,3);kSigma0(3,9)=kSigma0(9,3);kSigma0(3,10)=kSigma0(10,3);kSigma0(3,11)=kSigma0(11,3);kSigma0(3,12)=kSigma0(12,3);kSigma0(4,5)=kSigma0(5,4);
kSigma0(4,6)=kSigma0(6,4);kSigma0(4,7)=kSigma0(7,4);kSigma0(4,8)=kSigma0(8,4);kSigma0(4,9)=kSigma0(9,4);kSigma0(4,10)=kSigma0(10,4);kSigma0(4,11)=kSigma0(11,4);
kSigma0(4,12)=kSigma0(12,4);kSigma0(5,6)=kSigma0(6,5);kSigma0(5,7)=kSigma0(7,5);kSigma0(5,8)=kSigma0(8,5);kSigma0(5,9)=kSigma0(9,5);kSigma0(5,10)=kSigma0(10,5);
kSigma0(5,11)=kSigma0(11,5);kSigma0(5,12)=kSigma0(12,5);kSigma0(6,7)=kSigma0(7,6);kSigma0(6,8)=kSigma0(8,6);kSigma0(6,9)=kSigma0(9,6);kSigma0(6,10)=kSigma0(10,6);
kSigma0(6,11)=kSigma0(11,6);kSigma0(6,12)=kSigma0(12,6);kSigma0(7,8)=kSigma0(8,7);kSigma0(7,9)=kSigma0(9,7);kSigma0(7,10)=kSigma0(10,7);kSigma0(7,11)=kSigma0(11,7);
kSigma0(7,12)=kSigma0(12,7);kSigma0(8,9)=kSigma0(9,8);kSigma0(8,10)=kSigma0(10,8);kSigma0(8,11)=kSigma0(11,8);kSigma0(8,12)=kSigma0(12,8);kSigma0(9,10)=kSigma0(10,9);
kSigma0(9,11)=kSigma0(11,9);kSigma0(9,12)=kSigma0(12,9);kSigma0(10,11)=kSigma0(11,10);kSigma0(10,12)=kSigma0(12,10);kSigma0(11,12)=kSigma0(12,11);%12×12上三角对称矩阵赋值
kSigma0=N0*kSigma0;
%% ========================================================================
%应力刚度矩阵kSigma11
kSigma11=zeros(12);%初始化初应力刚度阵
kSigma11(2,2)=-6/(5*l);
kSigma11(4,4)=-Iz/(l*A);
kSigma11(6,6)=-2*l/15;
kSigma11(6,2)=-1/10;
kSigma11(11,6)=1/30;
kSigma11(8,2)=-kSigma11(2,2);kSigma11(8,6)=-kSigma11(6,2);kSigma11(10,4)=-kSigma11(4,4);kSigma11(12,2)=kSigma11(6,2);kSigma11(12,8)=-kSigma11(6,2);
kSigma11(2,8)=-kSigma11(2,2);kSigma11(6,8)=-kSigma11(6,2);kSigma11(4,10)=-kSigma11(4,4);kSigma11(2,12)=kSigma11(6,2);kSigma11(8,12)=-kSigma11(6,2);
kSigma11=(u_i-u_j)*E*A/l*kSigma11;
%% ========================================================================
%应力刚度矩阵kSigma21
kSigma21=zeros(12);%初始化初应力刚度阵
kSigma21(2,2)=-6/(5*l);kSigma21(4,4)=-Iy/(l*A);kSigma21(5,5)=-2*l/15;kSigma21(6,2)=-1/10;kSigma21(11,5)=1/30;%对角线元素
kSigma21(5,3)=1/10;kSigma21(9,3)=6/(5*l);kSigma21(10,4)=Iy/(l*A);kSigma21(11,3)=1/10;kSigma21(11,9)=-1/10;%上三角矩阵元素
kSigma21(3,5)=1/10;kSigma21(3,9)=6/(5*l);kSigma21(4,10)=Iy/(l*A);kSigma21(3,11)=1/10;kSigma21(9,11)=-1/10;%下三角矩阵元素
kSigma21=(u_i-u_j)*E*A/l*kSigma21;
%% ========================================================================
%局部到全局坐标变换                                                                    
k=k0+kL+kSigma0+kSigma11+kSigma21;%局部坐标系下单元刚度矩阵
ll=(x_j-x_i)/l;
mm=(y_j-y_i)/l;
nn=(z_j-z_i)/l;
DD=sqrt(ll^2+mm^2);
if n==1                           %判断是否和全局Z轴方向相同
lambda=[0,0,1;
        0,1,0;
        -1,0,0;];
elseif n==-1
lambda=[0,0,-1;
        0,1,0;
        1,0,0;];
else
lambda=[     ll,       mm,  mm;
          -mm/DD,    -ll/DD,   0;
       -ll*nn/DD, -mm*nn/DD,   DD;];
end
T=[lambda,zeros(3),zeros(3),zeros(3);
   zeros(3),lambda,zeros(3),zeros(3);
   zeros(3),zeros(3),lambda,zeros(3);
   zeros(3),zeros(3),zeros(3),lambda;]; 
k=T'*k*T;
end

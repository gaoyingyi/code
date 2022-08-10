%功能：
% 输入梁的材料参数、初始状态预张力，全局坐标系下的坐标、局部坐标系下的节点位移
% 输出全局坐标系下两节点十二自由度考虑剪切变形和应力作用下梁的刚度矩阵
function k=ESM_B31(E,mu,r,N0,nodeCoordinateMatrix,nodeDisplacementMatrix)
% E：弹性模量
% mu：泊松比
% nodeCoordinateMatrix：维度为2×3全局坐标系下节点坐标矩阵
% nodeDisplacementMatrix：维度为12×1的局部节点位移矩阵
%% ========================================================================
ui=nodeDisplacementMatrix(1,1);%#ok
vi=nodeDisplacementMatrix(2,1);%#ok
wi=nodeDisplacementMatrix(3,1);%#ok
thetaxi=nodeDisplacementMatrix(4,1);%#ok
thetayi=nodeDisplacementMatrix(5,1);%#ok
thetazi=nodeDisplacementMatrix(6,1);%#ok
uj=nodeDisplacementMatrix(7,1);%#ok
vj=nodeDisplacementMatrix(8,1);%#ok
wj=nodeDisplacementMatrix(9,1);%#ok
thetaxj=nodeDisplacementMatrix(10,1);%#ok
thetayj=nodeDisplacementMatrix(11,1);%#ok
thetazj=nodeDisplacementMatrix(12,1);%#ok

xi=nodeCoordinateMatrix(1,1);
yi=nodeCoordinateMatrix(1,2);
zi=nodeCoordinateMatrix(1,3);
xj=nodeCoordinateMatrix(2,1);
yj=nodeCoordinateMatrix(2,2);
zj=nodeCoordinateMatrix(2,3);
L=sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);%梁单元的长度
%% ========================================================================
%线性刚度阵预备参数
G=E/(2*(1+mu));%计算剪切模量
%--------------------------截面属性-----------------------------------------
D=2*r;%计算截面直径
A=pi*r^2;%计算截面面积
I=pi*D^4/64;%计算圆截面惯性矩,非线性矩阵中用的是这个，对于圆截面来说yz都一样
Iy=I;%计算圆截面y惯性矩
Iz=I;%计算圆截面z惯性矩
J=pi*D^4/32;%计算圆截面极惯性矩
Ay=A/10*9;%y方向有效抗剪面积
Az=A/10*9;%Z方向有效抗剪面积
phaY=12*E*Iz/(G*Ay*L^2);%沿单元坐标y方向的横向力剪切影响系数
phaZ=12*E*Iy/(G*Az*L^2);%沿单元坐标z方向的横向力剪切影响系数

%% ========================================================================
%线性刚度矩阵
k0=zeros(12);%初始化线性刚度阵
k0(1,1)=A*E/L;
k0(2,2)=12*E*Iz/((1+phaY)*L^3);
k0(3,3)=12*E*Iy/((1+phaZ)*L^3);
k0(4,4)=G*J/L;
k0(5,5)=(4+phaZ)*E*Iy/((1+phaZ)*L);
k0(6,6)=(4+phaY)*E*Iz/((1+phaY)*L);
k0(5,3)=-6*E*Iy/((1+phaZ)*L^2);
k0(6,2)=6*E*Iz/((1+phaY)*L^2);
k0(11,5)=(2-phaZ)*E*Iy/((1+phaZ)*L);
k0(12,6)=(2-phaY)*E*Iz/((1+phaY)*L);%这里列出了矩阵中10个不同的值，其他的值用这个赋予。
%下三角矩阵的其他值
k0(7,7)=k0(1,1);k0(8,8)=k0(2,2);k0(9,9)=k0(3,3);k0(10,10)=k0(4,4);k0(11,11)=k0(5,5);k0(12,12)=k0(6,6);
k0(7,1)=-k0(1,1);k0(8,2)=-k0(2,2);k0(8,6)=-k0(6,2);k0(9,3)=-k0(3,3);k0(9,5)=-k0(5,3);k0(10,4)=-k0(4,4);k0(11,3)=k0(5,3);k0(11,9)=-k0(5,3);k0(12,2)=k0(6,2);k0(12,8)=-k0(6,2);
% 对称矩阵，给矩阵上三角赋值
k0(1,7)=-k0(1,1);k0(2,8)=-k0(2,2);k0(6,8)=-k0(6,2);k0(3,9)=-k0(3,3);k0(5,9)=-k0(5,3);k0(4,10)=-k0(4,4);k0(3,11)=k0(5,3);k0(9,11)=-k0(5,3);k0(2,12)=k0(6,2);k0(8,12)=-k0(6,2);
k0(3,5)=k0(5,3);k0(2,6)=k0(6,2);k0(5,11)=k0(11,5);k0(6,12)=k0(12,6);
%% ========================================================================
%f非线性刚度矩阵
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'KL');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama0');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama11');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama12');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama13');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama21');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama22');
load('E:\code\matlab\TGLibrary\myData\Kbeam.mat', 'Ksigama23');
knonlinear=KL+Ksigama0+Ksigama11+Ksigama12+Ksigama13+Ksigama21+Ksigama22+Ksigama23;
knonlinear=eval(knonlinear);%将符号变量转换成double类型，如果有未赋值的变量就会报错
k=k0+knonlinear;%局部坐标系下单元刚度矩阵
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
k=T'*k*T;
end

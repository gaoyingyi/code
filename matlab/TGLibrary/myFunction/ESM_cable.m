%输入：弹性模量，泊松比，横截面半径，节点坐标，节点位移
%节点坐标矩阵包含已知的当前增量步的节点坐标
%节点位移矩阵包括已知节点位移（位移边界）和未知节点位移（待求位移）
function K=ESM_cable(E,mu,r,N0,nodeCoordinateMatrix,nodeDisplacementMatrix)

A=pi*r^2;%#ok
 syms ui vi wi uj vj wj
ui=nodeDisplacementMatrix(1,1);
vi=nodeDisplacementMatrix(2,1);
wi=nodeDisplacementMatrix(3,1);
uj=nodeDisplacementMatrix(4,1);
vj=nodeDisplacementMatrix(5,1);
wj=nodeDisplacementMatrix(6,1);

xi=nodeCoordinateMatrix(1,1);
yi=nodeCoordinateMatrix(1,2);
zi=nodeCoordinateMatrix(1,3);
xj=nodeCoordinateMatrix(2,1);
yj=nodeCoordinateMatrix(2,2);
zj=nodeCoordinateMatrix(2,3);

L0=sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);%#ok 梁单元的长度
L=sqrt((xi-xj-ui+uj)^2+(yi-yj-vi+vj)^2+(zi-zj-wi+wj)^2);
u=(ui - uj - xi + xj)/L;%#ok
v=(vi - vj - yi + yj)/L;%#ok
w=(wi - wj - zi + zj)/L;%#ok

load('E:\code\matlab\TGLibrary\myData\Kcable.mat', 'K');
K=eval(K);%将符号变量转换成double类型，如果有未赋值的变量就会报错
end
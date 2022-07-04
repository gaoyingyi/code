%两节点六自由度索单元，考虑内张力T，太简单无需测试
%输入：T应该是一个函数
%输入：L是原长
%输出6×6刚度阵、质量阵
function k=cableElementStiffnessMatrix(deltaL,L)
%% 参数
global  ECable;
ECable=20000;%弹性模量
r=1;%圆截面半径
A=pi*r^2;%计算截面面积
T=ECable*A*deltaL/L;
T=0;
%% 刚度阵
kl=ECable*A/L.*[1,0,0,-1,0,0;
           0,0,0, 0,0,0;
           0,0,0, 0,0,0;
          -1,0,0, 1,0,0;
           0,0,0, 0,0,0;
           0,0,0, 0,0,0;];                                                 %弹性刚度矩阵,补零是因为除了轴向力，其他方向不受力
ks=T/L.*[eye(3),-eye(3);-eye(3),eye(3)];                                   %初应力刚度矩阵
k=kl+ks;                                                                   %含内张力的非线性空间索单元刚度阵
end
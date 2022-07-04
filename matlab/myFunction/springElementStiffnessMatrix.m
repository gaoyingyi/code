%两节点两自由度弹簧单元
%扩充yz方向
%输入: T弹簧内力
%输入： L弹簧原长
%输出6×6刚度阵，集中质量阵
function k=springElementStiffnessMatrix(x)
k0=1e-6;
k1=0;
k3=0;
k5=0;

nonlinearSpringStiffness=k0+k1*x+k3*x^3+k5*x^5;                            %定义非线性弹簧刚度

k=nonlinearSpringStiffness.*[1,0,0,-1,0,0;
                             0,0,0, 0,0,0;
                             0,0,0, 0,0,0;
                            -1,0,0, 1,0,0;
                             0,0,0, 0,0,0;
                             0,0,0, 0,0,0;];                               %弹簧刚度矩阵
end
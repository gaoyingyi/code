%两节点两自由度弹簧单元
%扩充yz方向
%输入: T弹簧内力
%输入： L弹簧原长
%输出6×6刚度阵，集中质量阵
function [k,m]=springElementStiffnessMatrix(x)


k0=80;
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

m=1e-6.*eye(6);                                                            %集中质量矩阵
end

% nonlinearSpringStiffness=k1+2*k2*(b-delta)/delta*(sqrt(a^2+b^2)/sqrt(b-delta)^2+a^2-1);
%Kspring=20+60*x^2;  
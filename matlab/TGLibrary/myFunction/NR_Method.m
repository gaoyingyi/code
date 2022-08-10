%newton-raphson法解非线性方程组
% 输入：Fai为非线性方程组，Fai(u)~K(u)*u-F=0，应该是经过处理的静力学问题方程组
% 输入：K为切线刚度矩阵，是Fai的雅可比矩阵
% 输入：u初值
function s=NR_Method(Fai,K,u,eps)
if(~exist('eps','var'))
    eps = 10e-6;  % 默认精度
end
x=u(1,1);
y=u(2,1);
z=u(3,1);
J=eval(K);
F=eval(Fai);
delta =-inv(J)*F;
while norm(delta)>=eps        %norm范数
    u=u+delta;
x=u(1,1);
y=u(2,1);
z=u(3,1);
J=eval(K);
F=eval(Fai);
delta =-inv(J)*F;
end
s=delta+u;
return

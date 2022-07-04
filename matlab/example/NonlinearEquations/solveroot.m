function [x,fval] = solveroot
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
fun = @root2d; %非线性方程组
rng default
x0 = rand(2,1);%定义初值
[x,fval] = fsolve(fun,x0,options);
end
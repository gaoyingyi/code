clc;clear;close all;
%% fsolve 用法
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
fun = @root2d;
rng default
x0 = rand(2,1);
[x,fval] = fsolve(fun,x0,options);
%%
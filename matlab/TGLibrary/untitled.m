clc;clear;close all

%示例代码
fimplicit(@(CapitalOmega,X)((0-(CapitalOmega)^2)*X+3/4*0.4783*X^3+5/8*10.4241*X^5)^2+(2*0.01*(CapitalOmega)*X)^2-(0.025*(CapitalOmega)^2)^2, [0,0.4,0,0.4]); 
%fimplicit( @(x,y)x.^2 + y.^2 - 1, [ -1, 1 ]); 此语句的运行结果和上面语句相同
axis equal
axis( [0,0.4,0,0.4] )

%((0-[CapitalOmega]^2)*X+3/4*0.4783*X^3+5/8*10.4241*X^5)^2+(2*0.01*(CapitalOmega)*X)^2-(0.025*(CapitalOmega)^2)^2
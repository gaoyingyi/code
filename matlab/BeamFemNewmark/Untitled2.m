% 悬臂梁的自由振动
rho =1;A = 1;E = 1;I = 1;L = 1;c1 = 0;c2 = 0;Ne = 100;
xx = 1/Ne:1/Ne:1;
[Ma,Ka] = Beam3(rho,A,E,I,L/Ne,Ne+1);
Ca = c1*Ma+c2*Ka;
F = @(t) 0;             % 外载荷为零，自由振动
% for forced vibration e.g. F = @(t) sin(2*t)
% ------------------------------------------------
Lambda = 1.8751/L;
omega = Lambda^2*sqrt(E*I/rho/A);
h1 = cosh(Lambda*xx)-cos(Lambda*xx)-(cos(Lambda*L)+cosh(Lambda*L))...
    /(sin(Lambda*L)+sinh(Lambda*L))*(sinh(Lambda*xx)-sin(Lambda*xx));
h2 = Lambda*(sinh(Lambda*xx)+sin(Lambda*xx))-(cos(Lambda*L)+cosh(Lambda*L))...
    /(sin(Lambda*L)+sinh(Lambda*L))*(cosh(Lambda*xx)-cos(Lambda*xx))*Lambda;
D0 = zeros(2*Ne,1);
D0(1:2:2*Ne-1) = h1;
D0(2:2:2*Ne) = h2;
V0 = zeros(2*Ne,1);
dt = 1e-3;T = 20;tt = 0:dt:T;
[D,V,A] = Newmark(Ma,Ca,Ka,F,D0,V0,dt,T);
%% 精确解
ExactSol = @(x,t) (cosh(Lambda*x)-cos(Lambda*x)-(cos(Lambda*L)+cosh(Lambda*L))...
    /(sin(Lambda*L)+sinh(Lambda*L))*(sinh(Lambda*x)-sin(Lambda*x))).*cos(omega*t);
%% 
dn = D(Ne-1,:);     % newmark中点位移
dV = V(Ne-1,:);     % newmark中点位移
dA = A(Ne-1,:);     % newmark中点位移
de = ExactSol(L/2,tt);  %精确解
hold on
plot(tt,dn,'r','linewidth',2);
plot(tt,dV,'b','linewidth',2);
plot(tt,dA,'g','linewidth',2);

%plot(tt,de,'r')
title('Response of Mid-Point')
xlabel('t')
ylabel('Displacement')
legend('Newmark','Exact')


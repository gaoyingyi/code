clear all;close all;
%% 参数
E = 7e10;
rho = 2700;
L = 0.67;
b = 0.062;
h = 0.0025;
A = b*h;
I = b*h^3/12;
nL = 100;
l = L/nL;
% 模态自由度
ndof = 2;
symsdof = ndof*nL;  
% 单元矩阵
Ke = E*I/l^3.*[12,6*l,-12,6*l;
             6*l,4*l^2,-6*l,2*l^2;
             -12,-6*l,12,-6*l;
             6*l,2*l^2,-6*l,4*l^2;];
Me = rho*A*l/420.*[156,22*l,54,-13*l;
                 22*l,4*l^2,13*l,-3*l*l;
                 54,13*l,156,-22*l;
                 -13*l,-3*l*l,-22*l,4*l*l;];
KK = zeros(2*(nL+1),2*(nL+1));
MM = zeros(2*(nL+1),2*(nL+1));
             
for i = 1:nL
    KK((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2) = KK((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2)+Ke;
    MM((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2) = MM((i-1)*2+1:(i+1)*2,(i-1)*2+1:(i+1)*2)+Me;
end

K0 = KK(3:end,3:end);
M0 = MM(3:end,3:end);

[V,D]=eig(K0,M0);

DD = diag(D);
freq = DD.^0.5./2./pi;

dt = 0.0001;
t = 0:dt:4;
nt = length(t);

beta = 0.5;%步长
alpha = 0.25*(0.5+beta)^2;
c0 = 1/alpha/(dt^2);
c1 = beta/alpha/dt;
c2 = 1/alpha/dt;
c3 = 1/2/alpha-1;
c4 = beta/alpha-1;
c5 = dt/2*(beta/alpha-2);
c6 = dt*(1-beta);
c7 = beta*dt;
K_hat = K0+c0.*M0;
inverse_K_hat = inv(K_hat);

% initial conditions
X0 = zeros(symsdof,nt);
V0 = zeros(symsdof,nt);
a0 = zeros(symsdof,nt);

Q0 = zeros(symsdof,nt);

inverse_V = inv(V);

%Newmark method
for i = 2:nt
    Q0(symsdof-1,i) = 0.8*sin(26.39*i*dt);
    Q = Q0(:,i);
    
    Xt = X0(:,i-1);
    Vt = V0(:,i-1);
    at = a0(:,i-1);
    %control force
    
    Q_hat = Q+M0*(c0.*Xt+c2.*Vt+c3.*at);
    Xtdt = K_hat\Q_hat;
    atdt = c0.*(Xtdt-Xt)-c2.*Vt-c3.*at;
    Vtdt = Vt+c6.*at+c7.*atdt;
    
    X0(:,i) = Xtdt;
    V0(:,i) = Vtdt;
    a0(:,i) = atdt;
end

figure(1);
% plot(t,X0(symsdof-1,:)-L*X0(symsdof,:))
plot(t,X0(symsdof-1,:))
s=plot(t,X0(symsdof-1,:));
set(s,'LineWidth', 1);
title('Newmark')
xlabel('时间 (s)')
ylabel('振幅 (m)')
figure(2);
plot(t,V0(symsdof-1,:))
s=plot(t,V0(symsdof-1,:));
set(s,'LineWidth', 1);
title('Newmark')
xlabel('时间 (s)')
ylabel('振幅 (m)')
figure(3);
plot(t,a0(symsdof-1,:))
s=plot(t,a0(symsdof-1,:));
set(s,'LineWidth', 1);
title('Newmark')
xlabel('时间 (s)')
ylabel('振幅 (m)')

figure(4);
plot(t,Q0(symsdof-1,:))
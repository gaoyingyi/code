function [D,V,A] = Newmark(M,C,K,F,D0,V0,dt,T,Beta,Gamma)
% Newmark method for linear time invariant system
% **************    M*Y''(t)+C*Y'(t)+K*Y(t)=F(t)    **********************
% D0  - initial displacement
% V0  - initial velocity
% dt  - time increment
% T   - final time 
if nargin < 9
    Beta = 1/4;
    Gamma = 1/2;
end %
%% 常量参数
c1 = 1/Beta/dt^2;
c2 = Gamma/Beta/dt;
c3 = 1/Beta/dt;
c4 = 1/2/Beta-1;
c5 = Gamma/Beta-1;
c6 = (Gamma/2/Beta-1)*dt;
c7 = (1-Gamma)*dt;
c8 = Gamma*dt;
% -------------------------------------------------------------------------
A0 = M\(F(0)-K*D0-C*V0);     %初始加速度
n = T/dt+1;
m = length(D0);
D = zeros(m,n);
V = zeros(m,n);
A = zeros(m,n);
D(:,1) = D0;
V(:,1) = V0;
A(:,1) = A0;
% -------------------------------------------------------------------------
Kbar = c1*M+c2*C+K;         %等效刚度
for i = 1:n-1
    Da = D(:,i);
    Va = V(:,i);
    Aa = A(:,i);
    Fbar = F(i*dt)+M*(c1*Da+c3*Va+c4*Aa)+C*(c2*Da+c5*Va+c6*Aa);
    D(:,i+1) = Kbar\Fbar;
    A(:,i+1) = c1*(D(:,i+1)-Da)-c3*Va-c4*Aa;
    V(:,i+1) = Va+c7*Aa+c8*A(:,i+1);
end
end
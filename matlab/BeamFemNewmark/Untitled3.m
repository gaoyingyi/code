rho =1;A = 1;E = 1;I = 1;L = 1;c1 = 0.1;c2 = 0.1;
Ne = 100;                                   % number of elements
xx = 1/Ne:1/Ne:1;
h1 = 3*xx.^2-7*xx.^3+11/2*xx.^4-3/2*xx.^5;  % initial displacements
h2 = 6*xx-21*xx.^2+22*xx.^3-15/2*xx.^4;     % initial rotations
D0 = zeros(2*Ne,1);
D0(1:2:2*Ne-1) = h1;
D0(2:2:2*Ne) = h2;
V0 = zeros(2*Ne,1);
[Ma,Ka] = Beam3(rho,A,E,I,L/Ne,Ne+1);
Ca = c1*Ma+c2*Ka;
dt = 1e-3;
T = 10;
F = @(t) ExternalForce(t,Ne); 
[D,V,A] = Newmark(Ma,Ca,Ka,F,D0,V0,dt,T);
% -------------------------------------------------------------------------
subplot(311)
plot(0:dt:T,D(Ne-1,:),'linewidth',2)
xlabel('t')
title('Displacement of mid-point')
subplot(312)
plot(0:dt:T,D(2*Ne-1,:),'linewidth',2)
xlabel('t')
title('Displacement of end-point')
subplot(313)
plot([0,xx],[0;D(1:2:2*Ne-1,end)],'linewidth',2)
xlabel('x')
title('Mode Shape at the end time')
% -------------------------------------------------------------------------

function f = ExternalForce(t,Ne)
f = zeros(2*Ne,1);
f(Ne-1) = 100*sin(2*pi*5*t);
end
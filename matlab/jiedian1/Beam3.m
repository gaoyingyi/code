function [Ma,Ka,omega,lambdaL] = Beam3(rho,A,E,I,Le,n)
% 欧拉梁
% 离散化
% n 节点数
% Le 每个单元的长度
% -------------------------------------------------------------------------
% element stiffness matrix
Me = rho*A*Le*[13/35,11*Le/210,9/70,-13*Le/420;
               11*Le/210,Le^2/105,13*Le/420,-Le^2/140;
               9/70,13*Le/420,13/35,-11*Le/210;
               -13*Le/420,-Le^2/140,-11*Le/210,Le^2/105];
Ke = E*I/Le^2*[12/Le,6,-12/Le,6;
               6,4*Le,-6,2*Le;
               -12/Le,-6,12/Le,-6;
               6,2*Le,-6,4*Le];
% -------------------------------------------------------------------------
% global stiffness matrix
Ma = zeros(2*n,2*n);
Ka = zeros(2*n,2*n);
for i = 1:2:2*n-3
    Ma(i:i+3,i:i+3) = Ma(i:i+3,i:i+3) + Me;
    Ka(i:i+3,i:i+3) = Ka(i:i+3,i:i+3) + Ke;
end
% -------------------------------------------------------------------------
% boundary conditions !
bcs = 'cantilever';
%bcs = 'simply-supported';
switch bcs
    case 'cantilever'
        % the left end is clamped !
        Ma(:,1:2) = [];Ma(1:2,:) = [];
        Ka(:,1:2) = [];Ka(1:2,:) = [];
    case 'simply-supported'
        % simply supported at two ends
        Ma(:,[1,end-1]) = [];Ma([1,end-1],:) = [];
        Ka(:,[1,end-1]) = [];Ka([1,end-1],:) = [];
end
% -------------------------------------------------------------------------
% natural frequency
[~,w2] = eig(Ka,Ma);
w2 = diag(w2);
omega = sqrt(w2);
lambdaL = sqrt(omega)*(rho*A/E/I)^(1/4);
end
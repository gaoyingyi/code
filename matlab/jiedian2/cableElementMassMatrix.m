function m=cableElementMassMatrix(L)
%% 参数
global rhoCable;
rhoCable=1.5e-9;%索密度
r=1;            %横截面半径
A=pi*r^2;       %横截面面积
%% 一致质量阵，已检查
m=rhoCable*A*L/6.*[2*eye(3),eye(3);eye(3),2*eye(3)];
%%集中质量阵
%m=rho*A*L/2.*eye(6);
end
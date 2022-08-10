function m=EMM_Cable(E,mu,rho,r,nodeCoordinateMatrix)
%% 参数
xi=nodeCoordinateMatrix(1,1);
yi=nodeCoordinateMatrix(1,2);
zi=nodeCoordinateMatrix(1,3);
xj=nodeCoordinateMatrix(2,1);
yj=nodeCoordinateMatrix(2,2);
zj=nodeCoordinateMatrix(2,3);
L=sqrt((xi-xj)^2+(yi-yj)^2+(zi-zj)^2);%梁单元的长度

A=pi*r^2;       %横截面面积

%% 一致质量阵
m=rho*A*L/6.*[2*eye(3),eye(3);eye(3),2*eye(3)];

%% ========================================================================
%局部到全局坐标变换                                                                    
ll=(xj-xi)/L;
mm=(yj-yi)/L;
nn=(zj-zi)/L;
DD=sqrt(ll^2+mm^2);
if nn==1                           %判断是否和全局Z轴方向相同
    lambda=[0,0,1;
            0,1,0;
            -1,0,0;];
    elseif nn==-1                   %判断是否和全局Z轴方向相反
    lambda=[0,0,-1;
            0,1,0;
            1,0,0;];
    else                            %一般情况
    lambda=[     ll,       mm,  mm;
              -mm/DD,    -ll/DD,   0;
           -ll*nn/DD, -mm*nn/DD,   DD;];
end
T=[lambda,zeros(3);
   zeros(3),lambda;]; 
m=T'*m*T;

end
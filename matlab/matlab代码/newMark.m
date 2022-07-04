%% 隐式有限元动力学求解，newmark法求解位移、速度、加速度场
%输入：去掉奇异之后的M质量矩阵、K刚度矩阵、C阻尼矩阵
%输入：deltaT%时间步长
%输入：T %总时间
%输出：每一个自由度的位移场U,速度场UD,加速度场UDD
%输出：每一个自由度的外力载荷随时间的变化
function [U,UD,UDD,sumQ]=newMark(M,K,C,F,deltaT,T)  
%% 参数
    delta=0.5;alpha=0.25;
    [m,~]=size(M);       %去掉奇异后有多少个自由度，行
    n = T/deltaT+1;      %时间离散的个数，计算输出有多少列，加1是计算着t=0时刻的个数，列
    c0=1/(alpha*deltaT^2);c1=delta/(alpha*deltaT);c2=1/(alpha*deltaT);c3=1/(2*alpha)-1;c4=delta/alpha-1;c5=deltaT/2*(delta/alpha-2);c6=deltaT*(1-delta);c7=delta*deltaT; 
    KHU=K+c0*M+c1*C;     %等效刚度阵，时不变
%    [L,D,P]=ldl(KHU);  %三角分解KHU 
 %% 测试自由振动悬臂梁用的初始条件,不用的话注释掉
%    Ne=100;L=1;
%    Lambda = 1.8751/L;
%    xx = 1/Ne:1/Ne:1;
%    h1 = cosh(Lambda*xx)-cos(Lambda*xx)-(cos(Lambda*L)+cosh(Lambda*L))...
%         /(sin(Lambda*L)+sinh(Lambda*L))*(sinh(Lambda*xx)-sin(Lambda*xx));
%    h2 = Lambda*(sinh(Lambda*xx)+sin(Lambda*xx))-(cos(Lambda*L)+cosh(Lambda*L))...
%         /(sin(Lambda*L)+sinh(Lambda*L))*(cosh(Lambda*xx)-cos(Lambda*xx))*Lambda;
%    Q0=F(0);                        sumQ = zeros(m,n);  sumQ(:,1)=Q0;   %初始外载荷
%    U0 = zeros(2*Ne,1);
%    U0(1:2:2*Ne-1) = h1;
%    U0(2:2:2*Ne) = h2;U = zeros(m,n);     U(:,1)=U0;   %初始位移场
%    UD0 = zeros(2*Ne,1); UD = zeros(m,n);   UD(:,1)=UD0;   %初始速度场
%    UDD0=lsqminnorm(M,Q0-C*UD0-K*U0);UDD = zeros(m,n); UDD(:,1)=UDD0;   %初始加速度场
 %% 测试施加集中力的悬臂梁用的初始条件,不用的话注释掉
%     Ne = 100;                             % number of elements
%     xx = 1/Ne:1/Ne:1;
%     h1 = 3*xx.^2-7*xx.^3+11/2*xx.^4-3/2*xx.^5;  % initial displacements
%     h2 = 6*xx-21*xx.^2+22*xx.^3-15/2*xx.^4;     % initial rotations
%     Q0=F(0);                        sumQ = zeros(m,n);  sumQ(:,1)=Q0;   %初始外载荷
%     U0 = zeros(2*Ne,1);
%     U0(1:2:2*Ne-1) = h1;
%     U0(2:2:2*Ne) = h2;   U = zeros(m,n);     U(:,1)=U0;   %初始位移场
%     UD0 = zeros(2*Ne,1); UD = zeros(m,n);   UD(:,1)=UD0;   %初始速度场 
%     UDD0=lsqminnorm(M,Q0-C*UD0-K*U0);UDD = zeros(m,n); UDD(:,1)=UDD0;   %初始加速度场
%%  最简单的初始条件 
    Q0=F(0);                      
    sumQ = zeros(m,n);  
    sumQ(:,1)=Q0;   %初始外载荷
    U0=zeros(m,1);                     U = zeros(m,n);     U(:,1)=U0;   %初始位移场
    UD0=zeros(m,1);                   UD = zeros(m,n);   UD(:,1)=UD0;   %初始速度场
    UDD0=lsqminnorm(M,Q0-C*UD0-K*U0);UDD = zeros(m,n); UDD(:,1)=UDD0;   %初始加速度场
    %% 参数计算量，不需要改变，自动计算
    for i = 1:n-1                                                       %按时间步长循环
        Q=F(i*deltaT);                                                  %i+deltaT是本循环的时间
        sumQ(:,i+1)=Q;
        QHU=Q+M*(c0*U(:,i)+c2*UD(:,i)+c3*UDD(:,i))+C*(c1*U(:,i)+c4*UD(:,i)+c5*UDD(:,i));%t+delta时刻的QHU
        U(:,i+1)=KHU\QHU;                                               %计算出下一步的位移，然后储存到U里面。
        UDD(:,i+1)=c0*(U(:,i+1)-U(:,i))-c2*UD(:,i)-c3*UDD(:,i);
        UD(:,i+1)=UD(:,i)+c6*UDD(:,i)+c7*UDD(:,i+1);
    end
end
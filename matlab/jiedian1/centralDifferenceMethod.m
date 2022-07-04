%有限元动力学
%显示积分，中心差分法
%代码未完成
function []=centralDifferenceMethod(M,K,C,Q)
%输入：
%输出：
%%
%% 初始条件
u=0;v=0;
a=inv(M)(Q-C*v-K*u);
%% 参数
t=2;
T=100;%时间域
deltaT=0.001;%时间间隔
c0=1/(deltaT)^2;
c1=1/(2*deltaT);
c2=2*c0;
c3=1/c2;
u_deltaT=u0-deltaT*v+c3*a;
M=c0.*M+c1.*C;
%% 计算
for i=1:n
    Q = Q-(C-c2*M)u-(c0*M-c1*C)*uT_deltaT
    aTDeltaT=inv(M)*Q;
    a=c0*(aT_deltaT-2*ut+aTDeltaT);
    v=c1*(-uT_delta+aTdeltaT);
end

end
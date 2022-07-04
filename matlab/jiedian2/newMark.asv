%% 隐式有限元动力学求解，newmark法求解位移、速度、加速度场
%输入：deltaT%时间步长
%输入：T %总时间
%输出：每一个自由度的位移场U,速度场UD,加速度场UDD
%输出：每一个自由度的外力载荷随时间的变化
function [U,UD,UDD,sumQ]=newMark(T,deltaT)
%% 参数
%--------------------------------------------------------------------------
global sumsystemDof;                                      %总的系统自由度数目
global sumBoundarynodeDof;                                %受约束的自由度数目
global massDampingCoefficient;                                 %质量阻尼系数
global stiffnessDampingCoefficient;                            %刚度阻尼系数
%--------------------------------------------------------------------------
delta=0.5;
alpha=0.25;
c0=1/(alpha*deltaT^2);
c1=delta/(alpha*deltaT);
c2=1/(alpha*deltaT);
c3=1/(2*alpha)-1;
c4=delta/alpha-1;
c5=deltaT/2*(delta/alpha-2);
c6=deltaT*(1-delta);
c7=delta*deltaT;
%--------------------------------------------------------------------------
m=sumsystemDof-sumBoundarynodeDof;              %根据未约束的自由度数计算行数
n=T/deltaT+1;                                            %根据时间步计算列数
%-----------------------------坐标变换矩阵----------------------------------
lambdaBeam=coordinateTransformation(0,0,0,400,0,0,0,300,0);
lambdacable=coordinateTransformation(320,60,0,0,300,0,0,0,0);
lambdaSpring=coordinateTransformation(400,0,0,320,60,0,320,0,0);
beamCTF=[lambdaBeam,zeros(3),zeros(3),zeros(3);
    zeros(3),lambdaBeam,zeros(3),zeros(3);
    zeros(3),zeros(3),lambdaBeam,zeros(3);
    zeros(3),zeros(3),zeros(3),lambdaBeam;];                   %梁的转换矩阵
cableCTF=[lambdacable,zeros(3);
    zeros(3),lambdacable;];                                    %索的转换矩阵
springCTF=[lambdaSpring,zeros(3);
    zeros(3),lambdaSpring;];                                 %弹簧的转换矩阵
%-------------------------------质量矩阵----------------- ------------------
meBeam=beamElementMassmatrix(400);
meCable=cableElementMassMatrix(400);
meSpring=springElementMassMatrix();
meBeam=beamCTF'*meBeam*beamCTF;                 %整体坐标系下梁的单元质量矩阵
meCable=cableCTF'*meCable*cableCTF;             %整体坐标系下索的单元质量矩阵
meSpring=springCTF'*meSpring*springCTF;       %整体坐标系下弹簧的单元质量矩阵
M=zeros(sumsystemDof);                                   %初始化系统质量矩阵
M(1:12,1:12)=M(1:12,1:12)+meBeam;
M(7:9,7:9)=M(7:9,7:9)+meSpring(1:3,1:3);
M(7:9,13:15)=M(7:9,13:15)+meSpring(1:3,4:6);
M(13:15,7:9)=M(13:15,7:9)+meSpring(4:6,1:3);
M(13:15,13:15)=M(13:15,13:15)+meSpring(4:6,4:6);
M(13:18,13:18)=M(13:18,13:18)+meCable;                  %给系统质量矩阵赋值
%% 初始化
%-------------------------------零时刻刚度矩阵------------------------------
keBeam=beamElementStiffnessMatrix(0,400);
keCable=cableElementStiffnessMatrix(0,400);
keSpring=springElementStiffnessMatrix(0);
keBeam=beamCTF'*keBeam*beamCTF;                     %整体坐标系下梁的刚度矩阵
keCable=cableCTF'*keCable*cableCTF;                 %整体坐标系下索的刚度矩阵
keSpring=springCTF'*keSpring*springCTF;             %整体坐标系下索的刚度矩阵
K=zeros(sumsystemDof);                                   %初始化系统刚度矩阵
K(1:12,1:12)=K(1:12,1:12)+keBeam;
K(7:9,7:9)=K(7:9,7:9)+keSpring(1:3,1:3);
K(7:9,13:15)=K(7:9,13:15)+keSpring(1:3,4:6);
K(13:15,7:9)=K(13:15,7:9)+keSpring(4:6,1:3);
K(13:15,13:15)=K(13:15,13:15)+keSpring(4:6,4:6);
K(13:18,13:18)=K(13:18,13:18)+keCable;
%-------------------------------阻尼矩阵----------------- ------------------
C = massDampingCoefficient*M+stiffnessDampingCoefficient*K;
%-------------------------------去除边界------------------------------------
Ka=K(7:15,7:15);%除去位移边界后的刚度矩阵
Ma=M(7:15,7:15);%除去位移边界后的质量矩阵
Ca=C(7:15,7:15);%除去位移边界后的阻尼矩阵
%--------------------------------------------------------------------------
sumQ = zeros(m,n);
F=@(t) force(t,m);
Q0=F(0);
sumQ(:,1)=Q0;                                %初始化外载荷，并存入零时刻外载荷
%--------------------------------------------------------------------------
U = zeros(m,n);
U0=zeros(m,1);
U(:,1)=U0;                                     %初始位移场，并存入零时刻位移场
%--------------------------------------------------------------------------
UD = zeros(m,n);
UD0=zeros(m,1);
UD(:,1)=UD0;                                   %初始速度场，并存入零时刻速度场
%--------------------------------------------------------------------------
UDD = zeros(m,n);
UDD0=Ma\(Q0-Ca*UD0-Ka*U0);
UDD(:,1)=UDD0;                                %初始加速度场,并存入零时刻速度场
%--------------------------------------------------------------------------
i=1;%用于计数
for t = deltaT:deltaT:T    %按时间步长循环,每次循环计算i+1，即t时刻对应时刻的值
    KHU=Ka+c0.*Ma+c1.*Ca;                                        %等效刚度阵
    Q=F(t-deltaT);                                   %i+deltaT是本循环的时间
    sumQ(:,i)=Q;
    if t==T
        Q=F(t);
        sumQ(:,i+1)=Q;
    end
    QHU=Q+Ma*(c0*U(:,i)+c2*UD(:,i)+c3*UDD(:,i))+Ca*(c1*U(:,i)+ ...
        c4*UD(:,i)+c5*UDD(:,i));%t+delta时刻的QHU   
    U(:,i+1)=KHU\QHU;                     %计算出下一步的位移，然后储存到U里面
    UDD(:,i+1)=c0*(U(:,i+1)-U(:,i))-c2*UD(:,i)-c3*UDD(:,i);
    UD(:,i+1)=UD(:,i)+c6*UDD(:,i)+c7*UDD(:,i+1);
%-----------------------计算i+1，即t时刻对应的刚度值-------------------------
    keBeam=beamElementStiffnessMatrix(norm(U(1:3,i+1)),400);
    keCable=cableElementStiffnessMatrix(norm(-U(7:9,i+1)),400);
    keSpring=springElementStiffnessMatrix(norm(U(7:9,i+1)-U(1:3,i+1)));
%-------------------------------坐标转换------------------------------------
    lambdaBeam=coordinateTransformation(0,0,0,400,0,0,0,300,0);
    lambdacable=coordinateTransformation(320,60,0,0,300,0,0,0,0);
    lambdaSpring=coordinateTransformation(400,0,0,320,60,0,320,0,0);
    beamCTF=[lambdaBeam,zeros(3),zeros(3),zeros(3);
        zeros(3),lambdaBeam,zeros(3),zeros(3);
        zeros(3),zeros(3),lambdaBeam,zeros(3);
        zeros(3),zeros(3),zeros(3),lambdaBeam;];               %梁的转换矩阵
    cableCTF=[lambdacable,zeros(3);
        zeros(3),lambdacable;];                                %索的转换矩阵
    springCTF=[lambdaSpring,zeros(3);
        zeros(3),lambdaSpring;];                             %弹簧的转换矩阵
    keBeam=beamCTF'*keBeam*beamCTF;                 %整体坐标系下梁的刚度矩阵
    keCable=cableCTF'*keCable*cableCTF;             %整体坐标系下索的刚度矩阵
    keSpring=springCTF'*keSpring*springCTF;       %整体坐标系下弹簧的刚度矩阵
%-----------------------------刚度矩阵组集----------------------------------
    K=zeros(sumsystemDof);%初始化系统刚度矩阵
    K(1:12,1:12)=K(1:12,1:12)+keBeam;
    K(7:9,7:9)=K(7:9,7:9)+keSpring(1:3,1:3);
    K(7:9,13:15)=K(7:9,13:15)+keSpring(1:3,4:6);
    K(13:15,7:9)=K(13:15,7:9)+keSpring(4:6,1:3);
    K(13:15,13:15)=K(13:15,13:15)+keSpring(4:6,4:6);
    K(13:18,13:18)=K(13:18,13:18)+keCable;
    Ka=K(7:15,7:15);                                %除去位移边界后的刚度矩阵
%-----------------------------阻尼矩阵组集----------------------------------
    C = massDampingCoefficient*M+stiffnessDampingCoefficient*K;
    Ca=C(7:15,7:15);                                %除去位移边界后的阻尼矩阵
    i=i+1;
end
end
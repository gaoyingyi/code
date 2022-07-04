%主要实现：利用正运动学求解工作空间。
%注意：(不能直接运行！)
%1、随机函数N决定lim_theta角度的离散数值范围，rand(N,1)表示新建N个0~1的随机数，1表示1列，如果1改成10则
%是10列数，每列都是N个。如在theta1中用限制为[-pi,pi]那么为了使得离散的lim_theta1也在[-pi,pi]之间
%，因此其值为-pi+2*pi*rand(N,1)。得出结论X=[a,b]时，则离散范围X=a+(b-a)*rand(N,1)；

%2、for循环的n表示绘制点图时点的个数；

%3、subs函数对syms变量赋值，求解末端{6}的坐标系原点[0;0;0]相对基座坐标{0}的位置，由于T06是4X4矩阵，
%因此在{6}原点坐标后补1，用[0;0;0;1]表示其原点位置。
if exist('T06') && exist('R')  
    point=15000;
    lim_theta1=R(1,1)+(R(1,2)-R(1,1))*rand(point,1);
    lim_theta2=R(2,1)+(R(2,2)-R(2,1))*rand(point,1);
    lim_theta3=R(3,1)+(R(3,2)-R(3,1))*rand(point,1);
    lim_theta4=R(4,1)+(R(4,2)-R(4,1))*rand(point,1);
    lim_theta5=R(5,1)+(R(5,2)-R(5,1))*rand(point,1);
    lim_theta6=R(6,1)+(R(6,2)-R(6,1))*rand(point,1);
 
    parfor n=1:1:point
        p=subs(T06,{theta1 theta2 theta3 theta4 theta5 theta6},{lim_theta1(n),lim_theta2(n), ...
            lim_theta3(n),lim_theta4(n),lim_theta5(n),lim_theta6(n)})*[0;0;0;1];
        q(n,:) = p;
    end
    q = double(q);
    qx = q((1:point),1);
    qy = q((1:point),2);
    qz = q((1:point),3);
    scatter3(qx,qy,qz,'r','.')%画工作空间图
else
    disp('---------ERROR----------')
end
% % %----------------------------------------------------------------------------------
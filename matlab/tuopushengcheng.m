%生成三向网格
clc;clear;
%% 参数
n=6;      %分段数目
D=15000;   %天线口径
F=5000;    %天线焦距
delta=0; %定义垂度
N=[];    %存储六分之一节点坐标
NSum=[]; %存储所有节点的坐标
alpha=60*pi/180;
N1=zeros(3,n+1);%存储第一排节点
trans=[D/(4*n),(sqrt(3)*D)/(4*n),0]';%生成内部节点平移变换矩阵
%% 生成第一排节点的坐标
for i=1:n
    N1(1,i+1)=(D/2)/n*i;
end
N=[N,N1];      
%% 生成内部节点的坐标,考虑悬垂效应
for i=1:n-1
    NIn=N1(:,2:n-i+1)+trans*i; %存储内部节点
    N=[N,NIn];
end      %生成内部节点的坐标，没有经过悬垂效应改变边界节点坐标
if delta~=0
    R1=(D^2+16*delta^2)/(32*delta);%悬垂弧的半径
    gama=2*acos((R1-delta)/R1);%悬垂弧度
    transs=[(sqrt(3)*D/4-delta+R1)*cos(alpha/2),(sqrt(3)*D/4-delta+R1)*sin(alpha/2),0]';%由全局坐标系到局部坐标系的平移向量
    rotZZ=[cos(-(pi/6+gama/2)),sin(-(pi/6+gama/2)),0;
      -sin(-(pi/6+gama/2)),cos(-(pi/6+gama/2)),0;
         0,          0,     1;];     %由全局坐标系到局部坐标系的旋转矩阵
    t=n+1;
    for i=1:n-1
        N3=[-R1*cos(gama*i/n),R1*sin(gama*i/n),0]';
        N3=rotZZ*N3+transs;
        t=t+(n-i);
        N(:,t)=N3;
    end
end%悬垂效应修改边界


%% 旋转六次生成所有节点的坐标
NSum=[NSum,N];
for i=1:5
   rotz=[cos(alpha*i),sin(alpha*i),0;
     -sin(alpha*i),cos(alpha*i),0;
         0,          0,     1;];
    N2=rotz*N(:,2:end);
    NSum=[NSum,N2];
end
%% 投影法生成Z坐标
[m,t]=size(NSum);
for i=1:t
    NSum(3,i)=(NSum(1,i)^2+NSum(2,i)^2)/(4*F);
end
%% 画点
hold on
  plot3(NSum(1,:),NSum(2,:),NSum(3,:),'o','MarkerFaceColor','#1FFFFF')
  [r,m]=size(NSum);
  for i=1:m
    text(NSum(1,i),NSum(2,i),NSum(3,i),num2str(i),'FontSize',8,'Color','blue')
  end
  axis equal
u=1;%用于标号计数的
for i=1:3*n*(n+1)+1
     for j=1:3*n*(n+1)+1
          if (NSum(1,i)-NSum(1,j))^2+(NSum(2,i)-NSum(2,j))^2<=(D/(2*n))^2+0.001&&(NSum(1,i)-NSum(1,j))^2+(NSum(2,i)-NSum(2,j))^2>=(D/(2*n))^2-0.001&&i<j
             x=[NSum(1,i),NSum(1,j)];
             y=[NSum(2,i),NSum(2,j)];
             z=[NSum(3,i),NSum(3,j)];
             plot3(x,y,z,'r');  %划线的
             text((max(NSum(1,i),NSum(1,j))+min(NSum(1,i),NSum(1,j)))/2,(max(NSum(2,i),NSum(2,j))+min(NSum(2,i),NSum(2,j)))/2,(max(NSum(3,i),NSum(3,j))+min(NSum(3,i),NSum(3,j)))/2,num2str(u),'FontSize',10,'Color','red')
             u=u+1;
          end
     end
 end


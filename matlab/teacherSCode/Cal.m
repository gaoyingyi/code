clc;clear all;
format long;
global X0;
global Y0;
global Z0;
global N;            %节点总个数
global M;            %单元总个数
global M_B;          %梁单元个数
global M_XY;         %一类接触单元个数
global M_Z;          %二类接触单元个数
%global FContact_XY;
%global FContact_Z;
load Coordinate.txt;
load Beam.txt;
load Beam_XY.txt;
load Beam_Z.txt;
load Force.txt;
X0=1e-3*Coordinate(:,2);         %！注意单位m
Y0=1e-3*Coordinate(:,3);
Z0=1e-3*Coordinate(:,4);
[N,N1]=size(Coordinate);
[M_B,M_B1]=size(Beam);
[M_XY,M_XY1]=size(Beam_XY);
[M_Z,M_Z1]=size(Beam_Z);
% M=M_B+M_XY+M_Z;
K=zeros(6*N,6*N);
%%
%将梁单元刚度矩阵对号入座到整体刚度矩阵
global E_beam;            %弹性模量
global v_beam;            %泊松比
global d_beam;            %横截面的直径
E_beam=350000;            
v_beam=0.3;
d_beam=0.01;             %单位m
A_beam=0.25*d_beam^2*pi;
for m=1:M_B
    i=Beam(m,2);
    j=Beam(m,3);
    k=Beam(m,4);
    [Kii,Kij,Kji,Kjj]=KBeam(i,j,k,E_beam,v_beam,d_beam);
    K(6*i-5:6*i,6*i-5:6*i)=K(6*i-5:6*i,6*i-5:6*i)+Kii;
    K(6*i-5:6*i,6*j-5:6*j)=K(6*i-5:6*i,6*j-5:6*j)+Kij;
    K(6*j-5:6*j,6*i-5:6*i)=K(6*j-5:6*j,6*i-5:6*i)+Kji;
    K(6*j-5:6*j,6*j-5:6*j)=K(6*j-5:6*j,6*j-5:6*j)+Kjj;
end
%将第一类接触单元刚度矩阵对号入座到整体刚度矩阵
% global E_bar;        %接触单元弹性模量
% global A_bar;        %接触单元横截面积
% E_bar=3.5e11;
d_contact=18e-6;         %接触面积的直径        
for m=1:M_XY
    i=Beam_XY(m,2);
    j=Beam_XY(m,3);
    k=Beam_XY(m,4);
    [Kii,Kij,Kji,Kjj]=KBeam(i,j,k,E_beam,v_beam,d_contact);
    K(6*i-5:6*i,6*i-5:6*i)=K(6*i-5:6*i,6*i-5:6*i)+Kii;
    K(6*i-5:6*i,6*j-5:6*j)=K(6*i-5:6*i,6*j-5:6*j)+Kij;
    K(6*j-5:6*j,6*i-5:6*i)=K(6*j-5:6*j,6*i-5:6*i)+Kji;
    K(6*j-5:6*j,6*j-5:6*j)=K(6*j-5:6*j,6*j-5:6*j)+Kjj;
end
%将第二类接触单元刚度矩阵对好入座到整体刚度举证
for m=1:M_Z
    i=Beam_Z(m,2);
    j=Beam_Z(m,3);
    k=Beam_Z(m,4);
    [Kii,Kij,Kji,Kjj]=KBeam(i,j,k,E_beam,v_beam,d_contact);
    K(6*i-5:6*i,6*i-5:6*i)=K(6*i-5:6*i,6*i-5:6*i)+Kii;
    K(6*i-5:6*i,6*j-5:6*j)=K(6*i-5:6*i,6*j-5:6*j)+Kij;
    K(6*j-5:6*j,6*i-5:6*i)=K(6*j-5:6*j,6*i-5:6*i)+Kji;
    K(6*j-5:6*j,6*j-5:6*j)=K(6*j-5:6*j,6*j-5:6*j)+Kjj;
end
%%
%消除整体刚度矩阵的奇异
[nf,nf1]=size(Force);
BounNode=Force(:,1);
for i=1:nf
    for j=1:6
        if(0==Force(i,j+1))
            K(:,BounNode(i)*6-6+j)=0;
            K(BounNode(i)*6-6+j,:)=0;
            K(BounNode(i)*6-6+j,BounNode*6-6+j)=1;
        end
    end
end 

%%
%组集外力列向量
%[nf,nf1]=size(Force);
F=zeros(6*N,1);
for n=1:nf
    FNode=Force(n,1);
    for nn=1:6
       F(6*FNode-6+nn)=Force(n,nn+1);
    end
end

%%
%计算位移，提取接触力
u=inv(K)*F;
%Guass法求解
%u=Guass(K,F);
%变形后节点位置
for n=1:N
    X(n,1)=X0(n,1)+u(6*n-5);
    Y(n,1)=Y0(n,1)+u(6*n-4);
    Z(n,1)=Z0(n,1)+u(6*n-3);
end
%提取接触力
fid=fopen('result.txt','wt');
fprintf(fid,'%s\n','第一类接触力（XY平面），正数受压负数受拉:');
fprintf(fid,'%s\n','单元编号  接触应变 接触应力Pa  接触力N');
for n=1:M_XY
    i=Beam_XY(n,2);
    j=Beam_XY(n,3);
    l0=sqrt([X0(j)-X0(i)]^2+[Y0(j)-Y0(i)]^2+[Z0(j)-Z0(i)]^2);
    l=sqrt([X(j)-X(i)]^2+[Y(j)-Y(i)]^2+[Z(j)-Z(i)]^2);
    AL0=[X0(j)-X0(i);Y0(j)-Y0(i);Z0(j)-Z0(i)];
    AL=[X(j)-X(i);Y(j)-Y(i);Z0(j)-Z0(i)];

    FContact_XY(n,1)=Beam_XY(n,1);
    FContact_XY(n,2)=abs(l0-l)/l0;
    FContact_XY(n,3)=abs(l0-l)/l0*E_beam;
    FContact_XY(n,4)=abs(l0-l)/l0*E_beam*A_beam;
%    CoslTol0=AL0'*AL/l0/l;
%    FContact_XY(n,2)=(l0-l*CoslTol0)/l0;
%    FContact_XY(n,3)=(l0-l*CoslTol0)/l0*E_bar;
%    FContact_XY(n,4)=(l0-l*CoslTol0)/l0*E_bar*A_bar;
    fprintf(fid,'%d%s%1.4e%s%1.4e%s%f\n',FContact_XY(n,1),'   ',FContact_XY(n,2),'   ',FContact_XY(n,3),'   ',FContact_XY(n,4));
end
fprintf(fid,'\n\n%s\n','第二类接触力（Z轴），正数受压负数受拉:');
fprintf(fid,'%s\n','单元编号 接触应变 接触应力Pa  接触力N');
for n=1:M_Z
    i=Beam_Z(n,2);
    j=Beam_Z(n,3);
    l0=sqrt([X0(j)-X0(i)]^2+[Y0(j)-Y0(i)]^2+[Z0(j)-Z0(i)]^2);
    l=sqrt([X(j)-X(i)]^2+[Y(j)-Y(i)]^2+[Z(j)-Z(i)]^2);
    %l0=sqrt([Z0(j)-Z0(i)]^2);
    %l=sqrt([Z(j)-Z(i)]^2);
    FContact_Z(n,1)=Beam_Z(n,1);
    FContact_Z(n,2)=abs(l0-l)/l0;
    FContact_Z(n,3)=abs(l0-l)/l0*E_beam;
    FContact_Z(n,4)=abs(l0-l)/l0*E_beam*A_beam;
    fprintf(fid,'%d%s%1.4e%s%1.4e%s%f\n',FContact_Z(n,1),'   ',FContact_Z(n,2),'   ',FContact_Z(n,3),'   ',FContact_Z(n,4));
end
fclose(fid);

%%
%结果处理
for i=1:N
   Displacement(i)=sqrt(u(6*i-5)^2+u(6*i-4)^2+u(6*i-3)^2);
end
str=sprintf('最大位移量为：%e ，最小位移量为：%e',max(Displacement),min(Displacement))
%%
%绘图
hold on;
%变形前的图形
for n=1:M_B
    a=[X0(Beam(n,2)),X0(Beam(n,3))];
    b=[Y0(Beam(n,2)),Y0(Beam(n,3))];
    c=[Z0(Beam(n,2)),Z0(Beam(n,3))];
    plot3(a,b,c);
end
for n=1:M_XY
    a=[X0(Beam_XY(n,2)),X0(Beam_XY(n,3))];
    b=[Y0(Beam_XY(n,2)),Y0(Beam_XY(n,3))];
    c=[Z0(Beam_XY(n,2)),Z0(Beam_XY(n,3))];
    plot3(a,b,c);
end
for n=1:M_Z
    a=[X0(Beam_Z(n,2)),X0(Beam_Z(n,3))];
    b=[Y0(Beam_Z(n,2)),Y0(Beam_Z(n,3))];
    c=[Z0(Beam_Z(n,2)),Z0(Beam_Z(n,3))];
    plot3(a,b,c);
end
%变形后的图形
for n=1:M_B
    a=[X(Beam(n,2)),X(Beam(n,3))];
    b=[Y(Beam(n,2)),Y(Beam(n,3))];
    c=[Z(Beam(n,2)),Z(Beam(n,3))];
    plot3(a,b,c,'r');
end
for n=1:M_XY
    a=[X(Beam_XY(n,2)),X(Beam_XY(n,3))];
    b=[Y(Beam_XY(n,2)),Y(Beam_XY(n,3))];
    c=[Z(Beam_XY(n,2)),Z(Beam_XY(n,3))];
    plot3(a,b,c,'r');
end
for n=1:M_Z
    a=[X(Beam_Z(n,2)),X(Beam_Z(n,3))];
    b=[Y(Beam_Z(n,2)),Y(Beam_Z(n,3))];
    c=[Z(Beam_Z(n,2)),Z(Beam_Z(n,3))];
    plot3(a,b,c,'r');
end
hold off;

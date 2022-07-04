clc;
clear;
format short;
%%
% GUI用户输入
global F0;                  %正弦激励振幅   F=F0*sin(2*pi*f)
global f;                   %振动频率
global DirTag;              %基础振动加载方向 0-X方向，1-Y方向，2—Z方向
global VTypeTag;            %  0-正弦振动，1-随机振动
F0=0.001;
DirTag=2;
f=1000;
VTypeTag=1;
gmax=1e-2;           %功率谱最大幅值
fa=5;                %功率谱频率转角处值
fb=10;
fc=20;
fd=30;
time=5;              %加载时间
%%
global X0;
global Y0;
global Z0;
global N;            %节点总个数
global M;            %单元总个数
global M_B;          %梁单元个数
global M_XY;         %一类接触单元个数
global M_Z;          %二类接触单元个数
global density_beam;        %梁单元密度 
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
K=zeros(6*N,6*N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%模态分析%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%将梁单元刚度矩阵对号入座到整体刚度矩阵
global E_beam;            %弹性模量
global v_beam;            %泊松比
global d_beam;            %横截面的直径
E_beam=3.5e11;            
v_beam=0.28;
d_beam=52e-6;             %单位m
density_beam=19.35e3;
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
global d_Contact;        %接触单元横截面积
d_Contact=18e-6;         %接触面积     
for m=1:M_XY
    i=Beam_XY(m,2);
    j=Beam_XY(m,3);
    k=Beam_XY(m,4);
    [Kii,Kij,Kji,Kjj]=KBeam(i,j,k,E_beam,v_beam,d_Contact);
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
    [Kii,Kij,Kji,Kjj]=KBeam(i,j,k,E_beam,v_beam,d_Contact);
    K(6*i-5:6*i,6*i-5:6*i)=K(6*i-5:6*i,6*i-5:6*i)+Kii;
    K(6*i-5:6*i,6*j-5:6*j)=K(6*i-5:6*i,6*j-5:6*j)+Kij;
    K(6*j-5:6*j,6*i-5:6*i)=K(6*j-5:6*j,6*i-5:6*i)+Kji;
    K(6*j-5:6*j,6*j-5:6*j)=K(6*j-5:6*j,6*j-5:6*j)+Kjj;
end
%边界约束条件及载荷处理
[nf,nf1]=size(Force);
Force=sortrows(Force);
BounNode=Force(:,1);
InterNum=(1:6*N)';
OutsideNum=[];
Kdr=[];
Krr=[];
%处理列
for i=nf:-1:1
    for j=6:-1:1  
         InterNum(BounNode(i)*6-6+j)=[];
         OutsideNum=[OutsideNum;BounNode(i)*6-6+j];
         Kdr(:,i*6-6+j)= K(:,BounNode(i)*6-6+j);
         K(:,BounNode(i)*6-6+j)=[];
    end
end 
%处理行
for i=nf:-1:1
    for j=6:-1:1
         Krr(i*6-6+j,:)=Kdr(BounNode(i)*6-6+j,:);
         K(BounNode(i)*6-6+j,:)=[];
         Kdr(BounNode(i)*6-6+j,:)=[];
    end
end 
%将梁单元质量矩阵对号入座到整体质量矩阵
M=zeros(6*N,6*N);
for m=1:M_B
    i=Beam(m,2);
    j=Beam(m,3);
    [Mii,Mij,Mji,Mjj]=MBeam(i,j,E_beam,v_beam,d_beam,density_beam);
    M(6*i-5:6*i,6*i-5:6*i)=M(6*i-5:6*i,6*i-5:6*i)+Mii;
    M(6*i-5:6*i,6*j-5:6*j)=M(6*i-5:6*i,6*j-5:6*j)+Mij;
    M(6*j-5:6*j,6*i-5:6*i)=M(6*j-5:6*j,6*i-5:6*i)+Mji;
    M(6*j-5:6*j,6*j-5:6*j)=M(6*j-5:6*j,6*j-5:6*j)+Mjj;
end
%将第一类接触单元质量矩阵对号入座到整体质量矩阵
for m=1:M_XY
    i=Beam_XY(m,2);
    j=Beam_XY(m,3);
    [Mii,Mij,Mji,Mjj]=MBeam(i,j,E_beam,v_beam,d_Contact,density_beam);
    M(6*i-5:6*i,6*i-5:6*i)=M(6*i-5:6*i,6*i-5:6*i)+Mii;
    M(6*i-5:6*i,6*j-5:6*j)=M(6*i-5:6*i,6*j-5:6*j)+Mij;
    M(6*j-5:6*j,6*i-5:6*i)=M(6*j-5:6*j,6*i-5:6*i)+Mji;
    M(6*j-5:6*j,6*j-5:6*j)=M(6*j-5:6*j,6*j-5:6*j)+Mjj;
end
%将第二类接触单元质量矩阵对好入座到整体质量矩阵
for m=1:M_Z
    i=Beam_Z(m,2);
    j=Beam_Z(m,3);
    [Mii,Mij,Mji,Mjj]=MBeam(i,j,E_beam,v_beam,d_Contact,density_beam);
    M(6*i-5:6*i,6*i-5:6*i)=M(6*i-5:6*i,6*i-5:6*i)+Mii;
    M(6*i-5:6*i,6*j-5:6*j)=M(6*i-5:6*i,6*j-5:6*j)+Mij;
    M(6*j-5:6*j,6*i-5:6*i)=M(6*j-5:6*j,6*i-5:6*i)+Mji;
    M(6*j-5:6*j,6*j-5:6*j)=M(6*j-5:6*j,6*j-5:6*j)+Mjj;
end
%消除整体质量矩阵的奇异
for i=nf:-1:1
    for j=6:-1:1
        M(:,BounNode(i)*6-6+j)=[];
        M(BounNode(i)*6-6+j,:)=[];
    end
end 
%求固有频率和正则阵型
[V1,D1]=eig(K,M);
DD=diag(D1);               
[nd,nd1]=size(DD);
ND=(1:nd)';
X=sortrows([DD,ND,V1']);
DD=X(:,1);                 %特征值,即固有频率的平方
ND=X(:,2);
V=(X(:,3:(nd+2)))';
dd=[];
for i=1:nd
    if(DD(i)>0 && DD(i)==conj(DD(i)))
        dd=[dd;DD(i)];
    end
end
dd=sqrt(dd)/2/pi;          %四周固定的结构固有频率
for i=1:nd
    V(:,i)=V(:,i)/sqrt(sum(V(:,i).^2));    %求正则阵型
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%模态分析结束%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%振动响应分析%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
switch(VTypeTag)
    case 0      %正弦振动
        %求等式右边力向量
        FBoun=zeros(6*nf,1);
        switch(DirTag)
            case 0
                for i=1:nf
                    FBoun(6*i-5)=F0;
                end
            case 1
                for i=1:nf
                    FBoun(6*i-4)=F0;
                end
            case 2
                for i=1:nf
                    FBoun(6*i-3)=F0;
                end
        end
        F=-Kdr*(Krr\FBoun);
        %主坐标系下力的振幅
        F=V'*F;
        %求主坐标系下的响应(阻尼比为0)
        Z=F./DD./abs(1-(2*pi*f./sqrt(DD)));
        %系统物理坐标系下的稳态响应
        X=V*Z;
        %节点位移u=U*sin(2*pi*t)
        U(InterNum,1)=X;
        U(OutsideNum,1)=Krr\FBoun;
   
    case 1     %随机振动
        %将功率谱转换成随机力[T0,F0]  
        [F0,T0]=psd(gmax,fa,fb,fc,fd,time);  %也可以改成高斯分布，用randn随机生成一组F
        FBoun=zeros(6*nf,1);
        [nt,nt1]=size(T0);
        U=zeros(nd,nt);
        Z=zeros(nd,nt);
        for i=1:nt
            switch(DirTag)
                case 0
                    for i=1:nf
                        FBoun(6*i-5)=F0(i);
                    end
                case 1
                    for i=1:nf
                        FBoun(6*i-4)=F0(i);
                    end
                case 2
                    for i=1:nf
                        FBoun(6*i-3)=F0(i);
                    end
            end
            F=-Kdr*(Krr\FBoun);
            %主坐标系下力的振幅
            F=V'*F;
            %求主坐标系下的响应(阻尼比为0)
            for j=1:i-1                              %杜哈美积分
                for k=1:nd
                    syms t;
                    h=1/dd(k)*sin(dd(k)*(T0(j)-t))*F0(j);
                    Z(k,i)=1/2*double(int(h,t,T0(j),T0(j+1)))+Z(k,i); %%不停的存储占用内存太大，最好直接绘图然后覆盖
                end
            end
            %系统物理坐标系下的稳态响应
            X(:,i)=V*Z;                 %%不停的存储占用内存太大，最好直接绘图然后覆盖
            %节点位移u=U*sin(2*pi*t)
            U(InterNum,i)=X(:,i);       %%不停的存储占用内存太大，最好直接绘图然后覆盖
            U(OutsideNum,i)=Krr\FBoun;  %%不停的存储占用内存太大，最好直接绘图然后覆盖
                
        end           
end

%%
%绘图u=U*sin(2*pi*t)
hold on
hold off
%
%%

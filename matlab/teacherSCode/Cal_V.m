clc;
clear;
format short;
%%
% GUI�û�����
global F0;                  %���Ҽ������   F=F0*sin(2*pi*f)
global f;                   %��Ƶ��
global DirTag;              %�����񶯼��ط��� 0-X����1-Y����2��Z����
global VTypeTag;            %  0-�����񶯣�1-�����
F0=0.001;
DirTag=2;
f=1000;
VTypeTag=1;
gmax=1e-2;           %����������ֵ
fa=5;                %������Ƶ��ת�Ǵ�ֵ
fb=10;
fc=20;
fd=30;
time=5;              %����ʱ��
%%
global X0;
global Y0;
global Z0;
global N;            %�ڵ��ܸ���
global M;            %��Ԫ�ܸ���
global M_B;          %����Ԫ����
global M_XY;         %һ��Ӵ���Ԫ����
global M_Z;          %����Ӵ���Ԫ����
global density_beam;        %����Ԫ�ܶ� 
load Coordinate.txt;
load Beam.txt;
load Beam_XY.txt;
load Beam_Z.txt;
load Force.txt;
X0=1e-3*Coordinate(:,2);         %��ע�ⵥλm
Y0=1e-3*Coordinate(:,3);
Z0=1e-3*Coordinate(:,4);
[N,N1]=size(Coordinate);
[M_B,M_B1]=size(Beam);
[M_XY,M_XY1]=size(Beam_XY);
[M_Z,M_Z1]=size(Beam_Z);
K=zeros(6*N,6*N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ģ̬����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%������Ԫ�նȾ���Ժ�����������նȾ���
global E_beam;            %����ģ��
global v_beam;            %���ɱ�
global d_beam;            %������ֱ��
E_beam=3.5e11;            
v_beam=0.28;
d_beam=52e-6;             %��λm
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
%����һ��Ӵ���Ԫ�նȾ���Ժ�����������նȾ���
global d_Contact;        %�Ӵ���Ԫ������
d_Contact=18e-6;         %�Ӵ����     
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
%���ڶ���Ӵ���Ԫ�նȾ���Ժ�����������նȾ�֤
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
%�߽�Լ���������غɴ���
[nf,nf1]=size(Force);
Force=sortrows(Force);
BounNode=Force(:,1);
InterNum=(1:6*N)';
OutsideNum=[];
Kdr=[];
Krr=[];
%������
for i=nf:-1:1
    for j=6:-1:1  
         InterNum(BounNode(i)*6-6+j)=[];
         OutsideNum=[OutsideNum;BounNode(i)*6-6+j];
         Kdr(:,i*6-6+j)= K(:,BounNode(i)*6-6+j);
         K(:,BounNode(i)*6-6+j)=[];
    end
end 
%������
for i=nf:-1:1
    for j=6:-1:1
         Krr(i*6-6+j,:)=Kdr(BounNode(i)*6-6+j,:);
         K(BounNode(i)*6-6+j,:)=[];
         Kdr(BounNode(i)*6-6+j,:)=[];
    end
end 
%������Ԫ��������Ժ�������������������
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
%����һ��Ӵ���Ԫ��������Ժ�������������������
for m=1:M_XY
    i=Beam_XY(m,2);
    j=Beam_XY(m,3);
    [Mii,Mij,Mji,Mjj]=MBeam(i,j,E_beam,v_beam,d_Contact,density_beam);
    M(6*i-5:6*i,6*i-5:6*i)=M(6*i-5:6*i,6*i-5:6*i)+Mii;
    M(6*i-5:6*i,6*j-5:6*j)=M(6*i-5:6*i,6*j-5:6*j)+Mij;
    M(6*j-5:6*j,6*i-5:6*i)=M(6*j-5:6*j,6*i-5:6*i)+Mji;
    M(6*j-5:6*j,6*j-5:6*j)=M(6*j-5:6*j,6*j-5:6*j)+Mjj;
end
%���ڶ���Ӵ���Ԫ��������Ժ�������������������
for m=1:M_Z
    i=Beam_Z(m,2);
    j=Beam_Z(m,3);
    [Mii,Mij,Mji,Mjj]=MBeam(i,j,E_beam,v_beam,d_Contact,density_beam);
    M(6*i-5:6*i,6*i-5:6*i)=M(6*i-5:6*i,6*i-5:6*i)+Mii;
    M(6*i-5:6*i,6*j-5:6*j)=M(6*i-5:6*i,6*j-5:6*j)+Mij;
    M(6*j-5:6*j,6*i-5:6*i)=M(6*j-5:6*j,6*i-5:6*i)+Mji;
    M(6*j-5:6*j,6*j-5:6*j)=M(6*j-5:6*j,6*j-5:6*j)+Mjj;
end
%���������������������
for i=nf:-1:1
    for j=6:-1:1
        M(:,BounNode(i)*6-6+j)=[];
        M(BounNode(i)*6-6+j,:)=[];
    end
end 
%�����Ƶ�ʺ���������
[V1,D1]=eig(K,M);
DD=diag(D1);               
[nd,nd1]=size(DD);
ND=(1:nd)';
X=sortrows([DD,ND,V1']);
DD=X(:,1);                 %����ֵ,������Ƶ�ʵ�ƽ��
ND=X(:,2);
V=(X(:,3:(nd+2)))';
dd=[];
for i=1:nd
    if(DD(i)>0 && DD(i)==conj(DD(i)))
        dd=[dd;DD(i)];
    end
end
dd=sqrt(dd)/2/pi;          %���̶ܹ��Ľṹ����Ƶ��
for i=1:nd
    V(:,i)=V(:,i)/sqrt(sum(V(:,i).^2));    %����������
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ģ̬��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%����Ӧ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
switch(VTypeTag)
    case 0      %������
        %���ʽ�ұ�������
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
        %������ϵ���������
        F=V'*F;
        %��������ϵ�µ���Ӧ(�����Ϊ0)
        Z=F./DD./abs(1-(2*pi*f./sqrt(DD)));
        %ϵͳ��������ϵ�µ���̬��Ӧ
        X=V*Z;
        %�ڵ�λ��u=U*sin(2*pi*t)
        U(InterNum,1)=X;
        U(OutsideNum,1)=Krr\FBoun;
   
    case 1     %�����
        %��������ת���������[T0,F0]  
        [F0,T0]=psd(gmax,fa,fb,fc,fd,time);  %Ҳ���Ըĳɸ�˹�ֲ�����randn�������һ��F
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
            %������ϵ���������
            F=V'*F;
            %��������ϵ�µ���Ӧ(�����Ϊ0)
            for j=1:i-1                              %�Ź�������
                for k=1:nd
                    syms t;
                    h=1/dd(k)*sin(dd(k)*(T0(j)-t))*F0(j);
                    Z(k,i)=1/2*double(int(h,t,T0(j),T0(j+1)))+Z(k,i); %%��ͣ�Ĵ洢ռ���ڴ�̫�����ֱ�ӻ�ͼȻ�󸲸�
                end
            end
            %ϵͳ��������ϵ�µ���̬��Ӧ
            X(:,i)=V*Z;                 %%��ͣ�Ĵ洢ռ���ڴ�̫�����ֱ�ӻ�ͼȻ�󸲸�
            %�ڵ�λ��u=U*sin(2*pi*t)
            U(InterNum,i)=X(:,i);       %%��ͣ�Ĵ洢ռ���ڴ�̫�����ֱ�ӻ�ͼȻ�󸲸�
            U(OutsideNum,i)=Krr\FBoun;  %%��ͣ�Ĵ洢ռ���ڴ�̫�����ֱ�ӻ�ͼȻ�󸲸�
                
        end           
end

%%
%��ͼu=U*sin(2*pi*t)
hold on
hold off
%
%%

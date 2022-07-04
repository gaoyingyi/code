function [Kii,Kij,Kji,Kjj]=KBeam(i,j,k,E,v,d)     %%梁单元刚度矩阵
global X0;
global Y0;
global Z0;
%单元坐标系下的单元刚度矩阵
Le=sqrt([X0(j)-X0(i)]^2+[Y0(j)-Y0(i)]^2+[Z0(j)-Z0(i)]^2);
Ix=pi/32*d^4;
Iy=pi/64*d^4;
Iz=pi/64*d^4;
G=E/2/(1+v);    
A=0.25*d*d*pi;

Ay=0.9*A;                                      %沿单元坐标y方向杆横截面的有效抗剪影响
Az=0.9*A;
pha_y=12*E*Iz/G/Ay/Le/Le;                      %沿单元坐标y方向的横向力剪切影响系数
pha_z=12*E*Iy/G/Az/Le/Le;                      %沿单元坐标z方向的横向力剪切影响系数

Kii=zeros(6);
Kii(1,1)=E*A/Le;
Kii(2,2)=12*E*Iz/(1+pha_y)/Le/Le/Le;
Kii(3,3)=12*E*Iy/(1+pha_z)/Le/Le/Le;
Kii(4,4)=G*Ix/Le;
Kii(5,5)=(4+pha_z)*E*Iy/(1+pha_z)/Le;
Kii(6,6)=(4+pha_y)*E*Iz/(1+pha_y)/Le;
Kii(6,2)=6*E*Iz/(1+pha_y)/Le/Le;
Kii(2,6)=Kii(6,2);
Kii(5,3)=-6*E*Iy/(1+pha_z)/Le/Le;
Kii(3,5)=Kii(5,3);

Kjj=zeros(6);
Kjj(1,1)=E*A/Le;
Kjj(2,2)=12*E*Iz/(1+pha_y)/Le/Le/Le;
Kjj(3,3)=12*E*Iy/(1+pha_z)/Le/Le/Le;
Kjj(4,4)=G*Ix/Le;
Kjj(5,5)=(4+pha_z)*E*Iy/(1+pha_z)/Le;
Kjj(6,6)=(4+pha_y)*E*Iz/(1+pha_y)/Le;
Kjj(6,2)=-6*E*Iz/(1+pha_y)/Le/Le;
Kjj(2,6)=Kjj(6,2);
Kjj(5,3)=6*E*Iy/(1+pha_z)/Le/Le;
Kjj(3,5)=Kjj(5,3);

Kji=zeros(6);
Kji(1,1)=-E*A/Le;
Kji(2,2)=-12*E*Iz/(1+pha_y)/Le/Le/Le;
Kji(3,3)=-12*E*Iy/(1+pha_z)/Le/Le/Le;
Kji(4,4)=-G*Ix/Le;
Kji(5,5)=(2-pha_z)*E*Iy/(1+pha_z)/Le;
Kji(6,6)=(2-pha_y)*E*Iz/(1+pha_y)/Le;
Kji(6,2)=6*E*Iz/(1+pha_y)/Le/Le;
Kji(2,6)=-Kji(6,2);
Kji(5,3)=-6*E*Iy/(1+pha_z)/Le/Le;
Kji(3,5)=-Kji(5,3);
Kij=(Kji)';

%坐标转换
lx=[X0(j)-X0(i)]/Le;
mx=[Y0(j)-Y0(i)]/Le;
nx=[Z0(j)-Z0(i)]/Le;

a=[Y0(j)-Y0(i)]*[Z0(k)-Z0(i)]-[Z0(j)-Z0(i)]*[Y0(k)-Y0(i)];
b=[Z0(j)-Z0(i)]*[X0(k)-X0(i)]-[X0(j)-X0(i)]*[Z0(k)-Z0(i)];
c=[X0(j)-X0(i)]*[Y0(k)-Y0(i)]-[Y0(j)-Y0(i)]*[X0(k)-X0(i)];
Det=sqrt(a^2+b^2+c^2);
lz=a/Det;
mz=b/Det;
nz=c/Det;
ly=mz*nx-nz*mx;
my=nz*lx-lz*nx;
ny=lz*mx-mz*lx;
T=[lx,mx,nx,0,0,0;ly,my,ny,0,0,0;lz,mz,nz,0,0,0;0,0,0,lx,mx,nx;0,0,0,ly,my,ny;0,0,0,lz,mz,nz];
Kii=T'*Kii*T;
Kij=T'*Kij*T;
Kji=T'*Kji*T;
Kjj=T'*Kjj*T;

end
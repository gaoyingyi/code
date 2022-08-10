%朱文,周水兴.空间梁单元显式切线刚度矩阵推导[J].公路交通技术,2008(05):40-49+53.
%通过此文件可得到梁单元的切线刚度矩阵
clc;clear;close all;
syms x y z L E G N0 A r theta
syms ui vi wi thetaxi thetayi thetazi uj vj wj thetaxj thetayj thetazj
N1=1-x/L;
N2=x/L;
N3=1-3*x^2/L^2+2*x^3/L^3;
N4=x-2*x^2/L+x^3/L^2;
N5=3*x^2/L^2-2*x^3/L^3;
N6=-x^2/L+x^3/L^2;
Nu=[N1,0,0,0,0,0,N2,0,0,0,0,0];
Nv=[0,N3,0,0,0,N4,0,N5,0,0,0,N6];
Nw=[0,0,N3,0,-N4,0,0,0,N5,0,-N6,0];
Ntheta=[0,0,0,N1,0,0,0,0,0,N2,0,0];
deltae=[ui vi wi thetaxi thetayi thetazi uj vj wj thetaxj thetayj thetazj].';
B0=[diff(Nu,'x',1);
    -y*diff(Nv,'x',2);
    -z*diff(Nw,'x',2);
    sqrt(z^2+y^2)*diff(Ntheta,'x',1);];
Bv=[(diff(Nv,'x',1)-y*diff(Ntheta,'x',1))*deltae*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1));
    0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0;];
Bw=[(diff(Nw,'x',1)+z*diff(Ntheta,'x',1))*deltae*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1));
    0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0;];
D=[E,0,0,0;
   0,E,0,0;
   0,0,E,0;
   0,0,0,G];
%弹性矩阵K0----------------------------------------------------------------- 
K0=B0.'*D*B0*r;
K0=subs(K0,str2sym('y'),str2sym('r*cos(theta)'));
K0=subs(K0,str2sym('z'),str2sym('r*sin(theta)'));
K0=int(K0,'x','0','L');
K0=int(K0,'r',str2sym('0'),str2sym('r'));
K0=int(K0,'theta',str2sym('0'),str2sym('2*pi'));
K0=subs(K0,str2sym('r^2*pi'),str2sym('A'));
K0=subs(K0,str2sym('A*r^2'),str2sym('4*I'));

K0=simplify(K0);
%大位移矩阵KL---------------------------------------------------------------
KL=B0.'*D*Bv+B0.'*D*Bw+Bv.'*D*B0+Bv.'*D*Bv+Bv.'*D*Bw+Bw.'*D*B0+Bw.'*D*Bv+Bw.'*D*Bw;
KL=KL*r;
KL=subs(KL,str2sym('y'),str2sym('r*cos(theta)'));
KL=subs(KL,str2sym('z'),str2sym('r*sin(theta)'));
KL=int(KL,'x','0','L');
KL=int(KL,'r',str2sym('0'),str2sym('r'));
KL=int(KL,'theta',str2sym('0'),str2sym('2*pi'));
KL=expand(KL);
KL=subs(KL,str2sym('r^2*pi'),str2sym('A'));
KL=subs(KL,str2sym('A*r^2'),str2sym('4*I'));

KL=simplify(KL);
%应力矩阵Ksigama11----------------------------------------------------------
Ksigama11=E*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*diff(Nu,'x',1)*deltae*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1));
Ksigama11=Ksigama11*r;
Ksigama11=subs(Ksigama11,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama11=subs(Ksigama11,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama11=int(Ksigama11,'x','0','L');
Ksigama11=int(Ksigama11,'r',str2sym('0'),str2sym('r'));
Ksigama11=int(Ksigama11,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama11=expand(Ksigama11);
Ksigama11=subs(Ksigama11,str2sym('r^2*pi'),str2sym('A'));
Ksigama11=subs(Ksigama11,str2sym('A*r^2'),str2sym('4*I'));

Ksigama11=simplify(Ksigama11);
%应力矩阵Ksigama12----------------------------------------------------------
Ksigama12=E/2*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*((diff(Nv,'x',1)-z ...
    *diff(Ntheta,'x',1))*deltae).'*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nv,'x',1)-z*diff(Ntheta,'x',1)));
Ksigama12=Ksigama12*r;
Ksigama12=subs(Ksigama12,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama12=subs(Ksigama12,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama12=int(Ksigama12,'x','0','L');
Ksigama12=int(Ksigama12,'r',str2sym('0'),str2sym('r'));
Ksigama12=int(Ksigama12,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama12=expand(Ksigama12);
Ksigama12=subs(Ksigama12,str2sym('r^2*pi'),str2sym('A'));
Ksigama12=subs(Ksigama12,str2sym('A*r^2'),str2sym('4*I'));

Ksigama12=simplify(Ksigama12);
%应力矩阵Ksigama13----------------------------------------------------------
Ksigama13=E/2*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*((diff(Nw,'x',1)+y ...
    *diff(Ntheta,'x',1))*deltae).'*((diff(Nw,'x',1)+y*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1)));
Ksigama13=Ksigama13*r;
Ksigama13=subs(Ksigama13,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama13=subs(Ksigama13,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama13=int(Ksigama13,'x','0','L');
Ksigama13=int(Ksigama13,'r',str2sym('0'),str2sym('r'));
Ksigama13=int(Ksigama13,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama13=expand(Ksigama13);
Ksigama13=subs(Ksigama13,str2sym('r^2*pi'),str2sym('A'));
Ksigama13=subs(Ksigama13,str2sym('A*r^2'),str2sym('4*I'));

Ksigama13=simplify(Ksigama13);
%应力矩阵Ksigama21----------------------------------------------------------
Ksigama21=E*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*diff(Nu,'x',1)*deltae*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1));
Ksigama21=Ksigama21*r;
Ksigama21=subs(Ksigama21,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama21=subs(Ksigama21,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama21=int(Ksigama21,'x','0','L');
Ksigama21=int(Ksigama21,'r',str2sym('0'),str2sym('r'));
Ksigama21=int(Ksigama21,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama21=expand(Ksigama21);
Ksigama21=subs(Ksigama21,str2sym('r^2*pi'),str2sym('A'));
Ksigama21=subs(Ksigama21,str2sym('A*r^2'),str2sym('4*I'));

Ksigama21=simplify(Ksigama21);
%应力矩阵Ksigama22----------------------------------------------------------
Ksigama22=E/2*((diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*((diff(Nv,'x',1) ...
    -z*diff(Ntheta,'x',1))*deltae).'*((diff(Nv,'x',1)-z*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)));
Ksigama22=Ksigama22*r;
Ksigama22=subs(Ksigama22,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama22=subs(Ksigama22,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama22=int(Ksigama22,'x','0','L');
Ksigama22=int(Ksigama22,'r',str2sym('0'),str2sym('r'));
Ksigama22=int(Ksigama22,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama22=expand(Ksigama22);
Ksigama22=subs(Ksigama22,str2sym('r^2*pi'),str2sym('A'));
Ksigama22=subs(Ksigama22,str2sym('A*r^2'),str2sym('4*I'));

Ksigama22=simplify(Ksigama22);
%应力矩阵Ksigama23----------------------------------------------------------
Ksigama23=E/2*((diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*((diff(Nw,'x',1)+y ...
    *diff(Ntheta,'x',1))*deltae).'*((diff(Nw,'x',1)+y*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)));
Ksigama23=Ksigama23*r;
Ksigama23=subs(Ksigama23,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama23=subs(Ksigama23,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama23=int(Ksigama23,'x','0','L');
Ksigama23=int(Ksigama23,'r',str2sym('0'),str2sym('r'));
Ksigama23=int(Ksigama23,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama23=expand(Ksigama23);
Ksigama23=subs(Ksigama23,str2sym('r^2*pi'),str2sym('A'));
Ksigama23=subs(Ksigama23,str2sym('A*r^2'),str2sym('4*I'));

Ksigama23=simplify(Ksigama23);
%初应力矩阵Ksigama0----------------------------------------------------------
Ksigama0=N0/A*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*(diff(Nv,'x',1)-y ...
    *diff(Ntheta,'x',1))+(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*(diff(Nw,'x',1) ...
    +z*diff(Ntheta,'x',1)));
Ksigama0=Ksigama0*r;
Ksigama0=subs(Ksigama0,str2sym('y'),str2sym('r*cos(theta)'));
Ksigama0=subs(Ksigama0,str2sym('z'),str2sym('r*sin(theta)'));
Ksigama0=int(Ksigama0,'x','0','L');
Ksigama0=int(Ksigama0,'r',str2sym('0'),str2sym('r'));
Ksigama0=int(Ksigama0,'theta',str2sym('0'),str2sym('2*pi'));
Ksigama0=expand(Ksigama0);
Ksigama0=subs(Ksigama0,str2sym('r^2*pi'),str2sym('A'));
Ksigama0=subs(Ksigama0,str2sym('A*r^2'),str2sym('4*I'));
Ksigama0=simplify(Ksigama0);
% 清除不需要的变量
save('E:\code\matlab\TGLibrary\myData\Kbeam.mat',"K0","KL",'Ksigama0',"Ksigama11",'Ksigama12','Ksigama13','Ksigama21',"Ksigama22","Ksigama23");
clearvars -except KL K0 Ksigama11 Ksigama12 Ksigama13 Ksigama21 Ksigama22 Ksigama23 Ksigama0;






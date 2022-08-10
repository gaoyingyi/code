%朱文,周水兴.空间梁单元显式切线刚度矩阵推导[J].公路交通技术,2008(05):40-49+53.
%通过此文件可得到梁单元的切线刚度矩阵
clc;clear;close all;
syms x y z L E G N0 A
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
K0=B0.'*D*B0;
K0=int(K0,'x','0','L');
K0=int(K0,'y',str2sym('-h/2'),str2sym('h/2'));
K0=int(K0,'z',str2sym('-W/2'),str2sym('W/2'));
K0=subs(K0,str2sym('h^3*W'),str2sym('12*Iz'));
K0=subs(K0,str2sym('W^3*h'),str2sym('12*Iy'));
K0=subs(K0,str2sym('W*h'),str2sym('A'));
K0=subs(K0,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
K0=simplify(K0);
%大位移矩阵KL---------------------------------------------------------------
KL=B0.'*D*Bv+B0.'*D*Bw+Bv.'*D*B0+Bv.'*D*Bv+Bv.'*D*Bw+Bw.'*D*B0+Bw.'*D*Bv+Bw.'*D*Bw;
KL=int(KL,'x','0','L');
KL=int(KL,'y',str2sym('-h/2'),str2sym('h/2'));
KL=int(KL,'z',str2sym('-W/2'),str2sym('W/2'));
KL=expand(KL);
KL=subs(KL,str2sym('h^3*W'),str2sym('12*Iz'));
KL=subs(KL,str2sym('W^3*h'),str2sym('12*Iy'));
KL=subs(KL,str2sym('h*W'),str2sym('A'));
KL=subs(KL,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
KL=simplify(KL);
%应力矩阵Ksigama11----------------------------------------------------------
Ksigama11=E*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*diff(Nu,'x',1)*deltae*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1));
Ksigama11=int(Ksigama11,'x','0','L');
Ksigama11=int(Ksigama11,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama11=int(Ksigama11,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama11=expand(Ksigama11);
Ksigama11=subs(Ksigama11,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama11=subs(Ksigama11,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama11=subs(Ksigama11,str2sym('h*W'),str2sym('A'));
Ksigama11=subs(Ksigama11,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama11=simplify(Ksigama11);
%应力矩阵Ksigama12----------------------------------------------------------
Ksigama12=E/2*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*((diff(Nv,'x',1)-z ...
    *diff(Ntheta,'x',1))*deltae).'*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nv,'x',1)-z*diff(Ntheta,'x',1)));
Ksigama12=int(Ksigama12,'x','0','L');
Ksigama12=int(Ksigama12,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama12=int(Ksigama12,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama12=expand(Ksigama12);
Ksigama12=subs(Ksigama12,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama12=subs(Ksigama12,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama12=subs(Ksigama12,str2sym('h*W'),str2sym('A'));
Ksigama12=subs(Ksigama12,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama12=simplify(Ksigama12);
%应力矩阵Ksigama13----------------------------------------------------------
Ksigama13=E/2*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*((diff(Nw,'x',1)+y ...
    *diff(Ntheta,'x',1))*deltae).'*((diff(Nw,'x',1)+y*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nv,'x',1)-y*diff(Ntheta,'x',1)));
Ksigama13=int(Ksigama13,'x','0','L');
Ksigama13=int(Ksigama13,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama13=int(Ksigama13,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama13=expand(Ksigama13);
Ksigama13=subs(Ksigama13,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama13=subs(Ksigama13,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama13=subs(Ksigama13,str2sym('h*W'),str2sym('A'));
Ksigama13=subs(Ksigama13,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama13=simplify(Ksigama13);
%应力矩阵Ksigama21----------------------------------------------------------
Ksigama21=E*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*diff(Nu,'x',1)*deltae*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1));
Ksigama21=int(Ksigama21,'x','0','L');
Ksigama21=int(Ksigama21,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama21=int(Ksigama21,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama21=expand(Ksigama21);
Ksigama21=subs(Ksigama21,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama21=subs(Ksigama21,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama21=subs(Ksigama21,str2sym('h*W'),str2sym('A'));
Ksigama21=subs(Ksigama21,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama21=simplify(Ksigama21);
%应力矩阵Ksigama22----------------------------------------------------------
Ksigama22=E/2*((diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*((diff(Nv,'x',1) ...
    -z*diff(Ntheta,'x',1))*deltae).'*((diff(Nv,'x',1)-z*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)));
Ksigama22=int(Ksigama22,'x','0','L');
Ksigama22=int(Ksigama22,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama22=int(Ksigama22,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama22=expand(Ksigama22);
Ksigama22=subs(Ksigama22,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama22=subs(Ksigama22,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama22=subs(Ksigama22,str2sym('h*W'),str2sym('A'));
Ksigama22=subs(Ksigama22,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama22=simplify(Ksigama22);
%应力矩阵Ksigama23----------------------------------------------------------
Ksigama23=E/2*((diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*((diff(Nw,'x',1)+y ...
    *diff(Ntheta,'x',1))*deltae).'*((diff(Nw,'x',1)+y*diff(Ntheta,'x',1)) ...
    *deltae)*(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)));
Ksigama23=int(Ksigama23,'x','0','L');
Ksigama23=int(Ksigama23,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama23=int(Ksigama23,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama23=expand(Ksigama23);
Ksigama23=subs(Ksigama23,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama23=subs(Ksigama23,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama23=subs(Ksigama23,str2sym('h*W'),str2sym('A'));
Ksigama23=subs(Ksigama23,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama23=simplify(Ksigama23);
%初应力矩阵Ksigama0----------------------------------------------------------
Ksigama0=N0/A*((diff(Nv,'x',1)-y*diff(Ntheta,'x',1)).'*(diff(Nv,'x',1)-y ...
    *diff(Ntheta,'x',1))+(diff(Nw,'x',1)+z*diff(Ntheta,'x',1)).'*(diff(Nw,'x',1) ...
    +z*diff(Ntheta,'x',1)));
Ksigama0=int(Ksigama0,'x','0','L');
Ksigama0=int(Ksigama0,'y',str2sym('-h/2'),str2sym('h/2'));
Ksigama0=int(Ksigama0,'z',str2sym('-W/2'),str2sym('W/2'));
Ksigama0=expand(Ksigama0);
Ksigama0=subs(Ksigama0,str2sym('h^3*W'),str2sym('12*Iz'));
Ksigama0=subs(Ksigama0,str2sym('W^3*h'),str2sym('12*Iy'));
Ksigama0=subs(Ksigama0,str2sym('h*W'),str2sym('A'));
Ksigama0=subs(Ksigama0,str2sym('A*(W^2+h^2)'),str2sym('12*J'));
Ksigama0=simplify(Ksigama0);
%清除不需要的变量
save('Kbeam.mat',"K0","KL",'Ksigama0',"Ksigama11",'Ksigama12','Ksigama13','Ksigama21',"Ksigama22","Ksigama23");
clearvars -except KL K0 Ksigama11 Ksigama12 Ksigama13 Ksigama21 Ksigama22 Ksigama23 Ksigama0;






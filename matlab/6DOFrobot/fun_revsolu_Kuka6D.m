function [S,Q] = fun_revsolu_Kuka6D(T)
%FUN_REVSOLU_KUKA6D 此处显示有关此函数的摘要
%   此处显示详细说明
r11=T(1,1);r12=T(1,2);r13=T(1,3);px=T(1,4);
r21=T(2,1);r22=T(2,2);r23=T(2,3);py=T(2,4);
r31=T(3,1);r32=T(3,2);r33=T(3,3);pz=T(3,4);

a1=500;a2=1300;a3=55;
d1=1045;d4=1025;

n=1;

theta1=atan(py/px);
c1=fun_revsolu_cos(theta1);s1=fun_revsolu_sin(theta1);

K=(px*px+py*py+(pz-d1)*(pz-d1)+a1*a1-a2*a2-a3*a3-d4*d4-2*a1*c1*px-2*a1*s1*py)/(2*a2);
KI = a3^2+d4^2-K^2;
if(KI <0 || abs(px)>2826 || abs(py)>2826 || pz >3700 || pz <-1045)
    S = 0;
    Q = [];
else
    for i=0:1
        theta3 = atan2(d4,a3)-atan2(K,(-1)^i*sqrt(abs(a3^2+d4^2-K^2)));
        c3=fun_revsolu_cos(theta3);s3=fun_revsolu_sin(theta3);

        ats23 = (a3-a2*s3)*(-c1*px-s1*py+a1)+(d4+a2*c3)*(-pz+d1);
        atc23 = (a3-a2*s3)*(-pz+d1)+(d4+a2*c3)*(c1*px+s1*py-a1);  
        theta23 = atan2(ats23,atc23);
        theta2 = theta23-theta3;
        c2=fun_revsolu_cos(theta2);s2=fun_revsolu_sin(theta2);
        s23=fun_revsolu_sin(theta2+theta3);c23=fun_revsolu_cos(theta2+theta3);

        ats4 = -r13*s1+r23*c1;
        atc4 = -r13*c1*s23-r23*s1*s23-r33*c23;
        for j = 0:1
            theta4 = atan2(ats4,atc4) + j*pi;
            c4=fun_revsolu_cos(theta4);s4=fun_revsolu_sin(theta4);

            ats5 = -r23*(c1*s4-c4*s1*s23)+r13*(s1*s4+c1*c4*s23)+c4*c23*r33;
            atc5 = c1*c23*r13-r33*s23+c23*r23*s1;
            for k =0:1
                theta5 =(-1)^k*atan2(ats5,atc5);
                c5=fun_revsolu_cos(theta5);s5=fun_revsolu_sin(theta5);

                ats6 = r22*(c5*(c1*s4-c4*s1*s23)+c23*s1*s5)-r12*(c5*(s1*s4+c1*c4*s23)-c1*c23*s5) - r32*(s5*s23 + c4*c5*c23);
                atc6 = r12*(c4*s1-c1*s4*s23)-r22*(c1*c4+s1*s4*s23)-c23*r32*s4;
                theta6 = atan2(ats6,atc6);

                Q(n,:) = fun_round_matrix([theta1 theta2 theta3 theta4 theta5 theta6]);
                n = n+1;
            end
        end
    end

    S = 0;

    for i = 1:8
        if(-3.14<Q(i,1) && Q(i,1)<3.14 && -2.267<Q(i,2) && Q(i,2)<0.349 && -1.744<Q(i,3) && Q(i,3)<2.512)
            if(-3.14<Q(i,4) && Q(i,4)<3.14 && -3.14<Q(i,5) && Q(i,5)<3.14 && -3.14<Q(i,6) && Q(i,6)<3.14)
                S = S+1;
            end 
        end
    end
end

end

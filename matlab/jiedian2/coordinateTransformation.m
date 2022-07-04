%% 局部坐标系到全局坐标系的转换矩阵
function [lambdas]=coordinateTransformation(x1,y1,z1,x2,y2,z2,x3,y3,z3)
L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
l=(x2-x1)/L;
m=(y2-y1)/L;
n=(z2-z1)/L;
A=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
B=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
C=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
G=sqrt(A^2+B^2+C^2);
H=sqrt((B*n-C*m)^2+(C*l-A*n)^2+(A*m-B*l)^2);
lambdas=[l,m,n;
        (B*n-C*m)/(H),(C*l-A*n)/(H),(A*m-B*l)/(H);
        A/G,B/G,C/G];
end
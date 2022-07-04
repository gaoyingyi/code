%newton-raphson法解非线性方程组，
% 输入：x为迭代初值，
% 输入：eps为允许误差
%需要给方程组fx1和其雅可比矩阵的dfx1
function s=newtonRaphsonMethod(x,eps)
F = fx1(x);%非线性方程组
J = -dfx1(x);    %非线性方程组导数,输出是雅可比矩阵
invJ = inv(J);
delta =invJ*F';
while norm(delta)>=eps%norm范数
    x=delta'+x;
    F=fx1(x);
    J=-dfx1(x);
    invJ=inv(J);
    delta=invJ*F';
end
s=delta'+x;
return

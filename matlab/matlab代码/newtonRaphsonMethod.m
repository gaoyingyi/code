%newton-raphson��������Է����飬
% ���룺xΪ������ֵ��
% ���룺epsΪ�������
%��Ҫ��������fx1�����ſɱȾ����dfx1
function s=newtonRaphsonMethod(x,eps)
F = fx1(x);%�����Է�����
J = -dfx1(x);    %�����Է����鵼��,������ſɱȾ���
invJ = inv(J);
delta =invJ*F';
while norm(delta)>=eps%norm����
    x=delta'+x;
    F=fx1(x);
    J=-dfx1(x);
    invJ=inv(J);
    delta=invJ*F';
end
s=delta'+x;
return

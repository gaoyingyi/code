%newton-raphson��������Է�����
% ���룺FaiΪ�����Է����飬Fai(u)~K(u)*u-F=0��Ӧ���Ǿ�������ľ���ѧ���ⷽ����
% ���룺KΪ���߸նȾ�����Fai���ſɱȾ���
% ���룺u��ֵ
function s=NR_Method(Fai,K,u,eps)
if(~exist('eps','var'))
    eps = 10e-6;  % Ĭ�Ͼ���
end
x=u(1,1);
y=u(2,1);
z=u(3,1);
J=eval(K);
F=eval(Fai);
delta =-inv(J)*F;
while norm(delta)>=eps        %norm����
    u=u+delta;
x=u(1,1);
y=u(2,1);
z=u(3,1);
J=eval(K);
F=eval(Fai);
delta =-inv(J)*F;
end
s=delta+u;
return

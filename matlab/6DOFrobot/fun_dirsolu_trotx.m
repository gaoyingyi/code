%x轴旋转矩阵
function Rx = fun_dirsolu_trotx(x,y)
if nargin > 1 && strcmp(y, 'd')
    x = x *pi/180;
end
Rx=[1 0 0 0;0 cos(x) -sin(x) 0;0 sin(x) cos(x) 0;0 0 0 1];
end
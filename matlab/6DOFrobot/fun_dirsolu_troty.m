%y轴旋转矩阵
function Ry = fun_dirsolu_troty(y,x)
if nargin > 1 && strcmp(x, 'd')
    y = y *pi/180;
end
Ry=[cos(y) 0 sin(y) 0;0 1 0 0;-sin(y) 0 cos(y) 0;0 0 0 1];
end
%z轴旋转矩阵
function Rz = fun_dirsolu_trotz(z,x)
if nargin > 1 && strcmp(x, 'd')
    z = z *pi/180;
end
Rz=[cos(z) -sin(z) 0 0;sin(z) cos(z) 0 0;0 0 1 0;0 0 0 1];
end
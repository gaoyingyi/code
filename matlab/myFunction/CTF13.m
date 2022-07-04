%% 一维到三维的坐标变换矩阵
%输入全局坐标系下的节点坐标
%输出刚度矩阵
function [ctf]=CTF13(x1,y1,z1,x2,y2,z2)
    L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    Cx=(x2-x1)/L;
    Cy=(y2-y1)/L;
    Cz=(z2-z1)/L;
    ctf=[Cx,Cy,Cz,0,0,0;
         0,0,0,Cx,Cy,Cz;];
end
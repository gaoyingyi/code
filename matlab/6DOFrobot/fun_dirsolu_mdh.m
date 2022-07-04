% % %主要实现：空间坐标变换，利用DH参数求解i-1到i绕自身按顺序的变化的矩阵。
% % %注意：2018版本MATLAB使用时，必须将pi定义为字符即包含syms pi的声明！！(有时候单独定义pi会报错，可以不加分号定义pi)

% % %例如仅仅只是定义了theta求变化矩阵时：
% % % syms theta;
% % % fun_dh(theta,0,0,pi/2)
% 结果为：(4967757600021511*sin(theta))/81129638414606681695789005144064,出错！！

% % %当把pi定义好后：
% % % syms theta;
% % % syms pi   %此时不加分号，防止报错。
% % % fun_dh(theta,0,0,pi/2)
% 结果为：[ cos(theta), -sin(theta),  0, 0]
% % % % % [          0,           0, -1, 0]
% % % % % [ sin(theta),  cos(theta),  0, 0]
% % % % % [          0,           0,  0, 1]
% % %----------------------------------------------------------------------------------

function T = fun_dirsolu_mdh(theta,d,a,alpha)
T = fun_dirsolu_trotx(alpha)*fun_dirsolu_transl(a,0,0)*fun_dirsolu_trotz(theta)*fun_dirsolu_transl(0,0,d);
end

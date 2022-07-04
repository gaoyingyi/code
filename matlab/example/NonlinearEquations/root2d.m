%返回F，如果要调用句柄：Fun=@root2d
function F = root2d(x)
F = zeros(2,1); % 分配返回数组
F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;
end
%% Gauss-Serdel迭代法的函数文件：  
function [y]=gauseidel(A,b,x0)  
ep=1.0e-6;
D=diag(diag(A));  
L=-tril(A,-1);  
U=-triu(A,1); 
B=(D-L)\U;  
f=(D-L)\b;  
y=B*x0+f;  
  
while norm(y-x0)>=ep  
    x0=y;
    y=B*x0+f;
end
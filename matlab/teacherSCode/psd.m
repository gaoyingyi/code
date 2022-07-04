%用户输入
function [F0,T0]=psd(gmax,fa,fb,fc,fd,time)
n=500;
%计算
det_f=(fd-fa)/n;
for i=1:n
    f(i)=fa+(i-1)*det_f+0.5*det_f;
    if(f(i)>fa && f(i)<=fb)
        G(i)=gmax/(fb-fa)*(f(i)-fa);
    elseif(f(i)<=fc)
        G(i)=gmax;
    else
        G(i)=gmax/(fc-fd)*(f(i)-fd);
    end
end
theta=random('Uniform',0,2*pi,1,i);
for j=1:n
F0(j,1)=sum(sqrt(2*det_f*G).*sin(2*pi*f*time/n*(j-1)+theta),2);
end
m=(1:j)';
T0=time/n*(m-1);
end
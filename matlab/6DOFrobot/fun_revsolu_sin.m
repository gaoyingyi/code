function s = fun_revsolu_sin(n)
d=sin(n);
if d<0.001
    s = 0;
else
    s=d;
end
end
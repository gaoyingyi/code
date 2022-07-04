function c = fun_revsolu_cos(n)
d=cos(n);
if d<0.001
    c = 0;
else
    c=d;
end
end
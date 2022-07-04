
function T = fun_round_matrix(A)
[a,b] = size(A);
for x = 1:a
    for y = 1:b
        T(x,y) = roundn(double(A(x,y)),-3);
    end
end
end
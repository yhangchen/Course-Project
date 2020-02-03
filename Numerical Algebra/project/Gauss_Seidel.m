function x = Gauss_Seidel(A,b,x0,v)
G_S = tril(A); x = x0;
for i = 1:v
    x = x + G_S\(b-A*x);
end
end
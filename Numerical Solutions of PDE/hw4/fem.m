% implement general FEM method
function u = fem(f,g,x, grid)
n = length(x);
A = spalloc(n, n, 3*n); F = zeros(n, 1);
A(1,1) = x(2) / 3 + 1 / (x(2) - x(1)) + 1 / x(1);
A(1,2) = (x(2) - x(1)) / 6 - 1 / (x(2) - x(1));
F(1) = x(1) * (f(0) + 2 * f(x(1))) / 6 + (x(2) - x(1)) * (f(x(2)) + 2 * f(x(1))) / 6;
for i = 2:n-1
    A(i,i) = (x(i+1) - x(i-1)) / 3 + 1 / (x(i+1) - x(i)) + 1 / (x(i) - x(i-1));
    A(i,i-1) = (x(i) - x(i-1)) / 6 - 1 / (x(i) - x(i-1));
    A(i,i+1) = (x(i+1) - x(i)) / 6 - 1 / (x(i+1) - x(i));
    F(i) = (x(i) - x(i-1)) * (f(x(i-1)) + 2 * f(x(i))) / 6 + (x(i+1) - x(i)) * (f(x(i+1)) + 2 * f(x(i))) / 6;
end
A(n, n-1) = (x(n) - x(n-1)) / 6 - 1 / (x(n) - x(n-1));
A(n,n) = 1 + (x(n) - x(n-1)) / 3 + 1 / (x(n) - x(n-1));
F(n) = (x(n) - x(n-1)) * f(x(n)) / 2 + g;
if strcmp(grid, 'uniform')
    u =  conjugate_gradient(A,F,1e-12,1e12);
else
    u = A\F;
end
end

function x = conjugate_gradient(A,b,epsilon,k_max)
x = b; k = 0; r = b - A*x; rho = r'*r;
while sqrt(rho) > norm(b)*epsilon & k < k_max
    k = k+1;
    if k == 1
        p = r;
    else
        beta = rho/rho1; p = r + beta*p;
    end
    w = A*p; alpha = rho/(p'*w); x = x + alpha*p;
    r = r - alpha*w; rho1 = rho; rho = r'*r;
end
end


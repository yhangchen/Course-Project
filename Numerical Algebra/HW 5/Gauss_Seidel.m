function [u, k, t] = Gauss_Seidel(A, b, epsilon, k_max)
U = -triu(A,1);
n = length(A);k = 0;
u = b;
init_err = norm(A*u-b);
err = init_err;
t0 = clock;
while err > epsilon*init_err & k < k_max
    u = U*u + b;
    for j = 1:n-1
        u(j)=u(j)/A(j,j);
        u(j+1:n)=u(j+1:n)-u(j)*A(j+1:n,j);
    end
    u(n)=u(n)/A(n,n);
    err = norm(A*u-b);
    k = k + 1;
end
t = etime(clock, t0);
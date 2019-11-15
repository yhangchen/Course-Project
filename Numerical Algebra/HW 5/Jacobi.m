function [u, k, t] = Jacobi(A, b, epsilon, k_max)
u = b; k = 0; D = diag(A);
init_err = norm(A*u-b);
err = init_err;
t0 = clock;
while err > epsilon*init_err & k < k_max
    u = u - (A * u - b)./D;
    err = norm(A*u-b);
    k = k + 1;
end
t = etime(clock, t0);
end
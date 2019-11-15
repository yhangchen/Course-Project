function [x,k,t] = conjugate_gradient(A,b,epsilon,k_max)
x = b; k = 0; r = b - A*x; rho = r'*r;
t0 = clock;
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
t = etime(clock, t0);
end

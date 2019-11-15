function [x,k,t] = preconditioned_conjugate_gradient(A,b,epsilon,k_max)
x = b; k = 0; r = b - A*x; normb = norm(b); d = diag(A);
t0 = clock;
while norm(r) > epsilon * normb & k < k_max
    z = r./d;
    k = k + 1;
    if k == 1
        p = z; rho = r'*z;
    else
        rho1 = rho; rho = r'*z; beta = rho/rho1; p = z + beta*p;
    end
    w = A*p; alpha = rho/(p'*w); x = x + alpha * p; r = r - alpha*w;
end
t = etime(clock, t0);
end

    
function [x, iter, out]=gl_cvx_mosek(x0, A, b, mu, opts)
iter = 0;
n = size(A, 2); l = size(b, 2);
cvx_solver mosek
cvx_begin quiet
variable x(n, l)
minimize(0.5* sum(sum(( A * x - b ).^2)) + mu * sum(norms(x,2,2)))
cvx_end
out = [];
out.val = 0.5* sum(sum(( A * x - b ).^2)) + mu * sum(norms(x,2,2));
end
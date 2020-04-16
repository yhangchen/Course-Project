function [X, out] = cvx_sol_1(S, mu)
% OUTPUT: X1: the solution
%         out: additional information.
[n,~] = size(S);
cvx_expert true
cvx_begin
    variable X(n,n) semidefinite
    maximize log_det(X) - trace(S*X) - mu * norm(vec(X), 1)
cvx_end
out.dualgap = n - trace(S*X)  - mu*norm(vec(X), 1);
out.value = log(det(X)) - trace(S*X) - mu * norm(vec(X), 1);
end
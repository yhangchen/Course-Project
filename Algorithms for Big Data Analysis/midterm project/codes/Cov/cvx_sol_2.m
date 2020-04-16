function [X, out] = cvx_sol_2(S, sigma)
% OUTPUT: X1: the solution
%         out: additional information.
[n,~] = size(S);
I = eye(n);
cvx_expert true
cvx_begin
    variable X(n,n) semidefinite
    minimize  norm(vec(X), 1) + sigma/2*sum(sum_square(vec(S*X-I)))
cvx_end
out.value = norm(vec(X), 1) + sigma/2*norm(S*X-I, 'fro')^2;
out.dualgap = norm(vec(X),1) + sigma * norm(S*X, 'fro')^2 - sigma * trace(S*X);
end
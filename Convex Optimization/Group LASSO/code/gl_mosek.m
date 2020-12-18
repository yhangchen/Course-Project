function [x,iter, out]=gl_mosek(x0, A, b, mu, opts)
iter=0;
clear prob;
[r, res] = mosekopt('symbcon');
[m, n] = size(A);
[~, l] = size(b);
A0 = sparse(A);
for i = 1:l-1; A0 = blkdiag(A0,A); end

% Linear constraints.
prob.a   = sparse(1,n*l+n+2);
prob.c = [zeros(n*l,1); 0; 1; mu*ones(n,1)];

% Affine conic constraint
tmp = sparse(2,n*l+n+2);
tmp(1,n*l+1) = 2;
tmp(2,n*l+2) = 1;
probF = sparse([zeros(1,n*l+1), 1, zeros(1, n);
    tmp]);
tmp = sparse([zeros(1,n*l), 1, zeros(1, 1+n);
    A0, zeros(m*l, n+2)]);
probF = [probF; tmp];
for i = 1:n
    tmp = spalloc(n*l+1,n*l+2+n,n+2);
    tmp(1,end-n+i) = 1;
    for j = 1:l
        tmp(i+1+(j-1)*n,i+(j-1)*n) = 1;
    end
    probF = [probF; tmp];
end
prob.f = probF;
prob.g = [0; 0; -4; 0; -vec(b); zeros((n*l+1)*n,1)]; % we replace t+2 with t
probCones = [res.symbcon.MSK_CT_QUAD 3 res.symbcon.MSK_CT_QUAD m*l+1];
for i = 1:n
    probCones = [probCones res.symbcon.MSK_CT_QUAD n*l+1];
end
prob.cones = probCones;
% Solve
[r, res] = mosekopt('minimize echo(0)', prob);
x = reshape(res.sol.itr.xx(1:n*l),n,l);
out.val = 0.5* sum(sum(( A * x - b ).^2)) + mu * sum(norms(x,2,2));
end
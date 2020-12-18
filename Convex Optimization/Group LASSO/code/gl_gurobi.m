function [x, iter, out]=gl_gurobi(x0, A, b, mu, opts)
iter=0;
[m, n] = size(A);
[~, l] = size(b);
A0 = sparse(A);
for i = 1:l-1; A0 = blkdiag(A0,A); end
model.obj = [zeros((m+n)*l,1); 0; 1; mu*ones(n,1)];
model.modelsense = 'min';
% transoform SOCP to QCQP.

% Add constraint: positive
model.A = [-speye(m*l,m*l), A0, sparse(m*l,n+2)];
model.rhs = vec(b);
model.sense = '=';
% specify each lower bound.
model.lb = [-inf*ones((m+n)*l,1); zeros(n+2,1)];
% Add second-order cone
A1 = sparse(2,(m+n)*l+n+2);
A1(1,(m+n)*l+1) = 2;
A1(2,(m+n)*l+2) = 1;
b1 = [0;-2];
c1 = [zeros(1,(m+n)*l+1), 1, zeros(1, n)];
d1 = 2;
model.quadcon(1).Qc = sparse(A1'*A1 - c1'*c1);
model.quadcon(1).q  = 2*(A1'*b1-c1'*d1);
model.quadcon(1).rhs = 0;

A2 = [speye(m*l,m*l), sparse(m*l,n*l+n+2)];
c2 = [zeros(1,(m+n)*l), 1, zeros(1, n+1)];
model.quadcon(2).Qc = sparse(A2'*A2 - c2'*c2);
model.quadcon(2).q = sparse((m+n)*l+n+2,1);
model.quadcon(2).rhs = 0;

for i = 1:n
    tmp = spalloc(n*l,(m+n)*l+2+n,l);
    for j = 1:l
        tmp(i+(j-1)*n,i+(j-1)*n+m*l) = 1;
    end
    A2 = tmp;
    c2 = sparse(1,(m+n)*l+n+2);
    c2(end-n+i) = 1;
    model.quadcon(i+2).Qc = sparse(A2'*A2 - c2'*c2);
    model.quadcon(i+2).q  = zeros((m+n)*l+2+n, 1);
    model.quadcon(i+2).rhs = 0;
end

result = gurobi(model);
x = reshape(result.x(m*l+1:(m+n)*l),n,l);
out.val = result.objval;


end
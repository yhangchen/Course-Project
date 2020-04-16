function [x,out]=l1l1_gurobi(x0, A, b, mu, opts)
	[m, n] = size(A);
    f1 = ones(2*n,1); f2 = ones(2*m,1);
    f = [mu*f1;f2];
    A1 = [A -A -eye(m) eye(m)];
    l = zeros(2*n+2*m,1);
	model = [];
	model.obj = f;
	model.A = sparse(A1);
	model.sense = '=';
	model.rhs = b;
    model.lb = l;
	results = gurobi(model);
	xx = results.x;
	x = xx(1:n)-xx(n+1:2*n);
	out = [];
	out.results = results;
	out.val = norm(A*x-b,1)+mu*norm(x,1);
end
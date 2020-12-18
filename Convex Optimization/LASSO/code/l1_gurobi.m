function [x,out]=l1_gurobi(x0, A, b, mu, opts)
	[m, n] = size(A);
	ATA = A'*A;
    ATb = A'*b;
	model = [];
	model.Q = 0.5*sparse([ATA,-ATA;-ATA,ATA]);
	model.obj = [-ATb;ATb]+mu*ones(2*n,1);
	model.A = speye(2*n);
    model.rhs = zeros(2*n,1);
	model.sense = '>';
	params.BarConvTol = 1e-12;
	info = gurobi(model, params);
	xx = info.x;
	x = xx(1:n)-xx(n+1:2*n);
	out = [];
	out.info = info;
	out.val = 0.5*norm(A*x-b)^2+mu*norm(x,1);
end
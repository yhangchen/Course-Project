function [x,out]=l1_mosek(x0, A, b, mu, opts)
	[~, n] = size(A);
	ATA = A'*A;
    ATb = A'*b;
    H = [ATA,-ATA;-ATA,ATA];
	f = [-ATb;ATb]+mu*ones(2*n,1);
	a = ones(1,2*n);
	[xx] = quadprog(H,f,[],[],[],[],sparse(2*n,1));
	x = xx(1:n)-xx(n+1:2*n);
	out = [];
	out.val = 0.5*norm(A*x-b)^2+mu*norm(x,1);
end
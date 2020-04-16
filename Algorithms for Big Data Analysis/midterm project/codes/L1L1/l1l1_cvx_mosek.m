function [x, out] = l1l1_cvx_mosek(x0, A, b, mu, opts)
    [~,n] = size(A);
	tic;
    cvx_solver mosek
	cvx_begin quiet
	    variable x(n)
	    minimize( norm( A * x - b, 1 ) + mu * norm( x, 1 ) )
	cvx_end
	time = toc;
	out = [];
	out.time = time;
	out.val = norm( A * x - b, 1 ) + mu * norm( x, 1 );
end 
function test1

	% min ||Ax-b||_1 + mu*||x||_1

	% generate data
	n = 1024;
	m = 512;

	% set the random seed
	rng(1234);

	A = randn(m,n);
	u = sprandn(n,1,0.1);
	b = A*u;
	mu = 1e-2;

	x0 = rand(n,1);

    errfun = @(x1, x2) norm(x1-x2,1)/(1+norm(x1,1));
    resfun = @(x) norm(A*x-b,1);
    nrm1fun = @(x) norm(x,1);
    
	% cvx calling mosek
	opts1 = []; %modify options
	tic; 
	[x1, out1] = l1l1_cvx_mosek(x0, A, b, mu, opts1);
	t1 = toc;

	% cvx calling gurobi
	opts2 = []; %modify options
	tic; 
	[x2, out2] = l1l1_cvx_gurobi(x0, A, b, mu, opts2);
	t2 = toc;
	% call mosek directly
	opts3 = []; %modify options
	tic; 
	[x3, out3] = l1l1_mosek(x0, A, b, mu, opts3);
	t3 = toc;

	% call gurobi directly
	opts4 = []; %modify options
	tic; 
	[x4, out4] = l1l1_gurobi(x0, A, b, mu, opts4);
	t4 = toc;
	% print comparison results with cvx-call-mosek
    fprintf('cvx-call-mosek: res: %3.2e, nrm1: %3.2e, cpu: %5.2f\n', resfun(x1), nrm1fun(x1), t1);
    fprintf('call-mosek: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x3), nrm1fun(x3), t3, errfun(x1, x3));
    fprintf('cvx-call-gurobi: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x2), nrm1fun(x2), t2, errfun(x1, x2));
    fprintf('call-gurobi: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x4), nrm1fun(x4), t4, errfun(x1, x4));
    
    fig = figure;
    hold on
    subplot(5,1,1)
    plot(u)
    legend('solution')
    subplot(5,1,2)
    plot(x1)
    legend('cvx-mosek')
    subplot(5,1,3)
    plot(x2)
    legend('cvx-gurobi')
    subplot(5,1,4)
    plot(x3)
    legend('mosek')
    subplot(5,1,5)
    plot(x4)
    legend('gurobi');
    hold off
end
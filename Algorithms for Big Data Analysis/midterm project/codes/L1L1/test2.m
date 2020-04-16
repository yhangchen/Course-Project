function test2

	% min ||Ax-b||_1 + mu*||x||_1

	% generate data
	n = 1024;
	m = 512;
    
	% set the random seed
    rng(1234)
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
    fprintf('cvx_mosek solved\n')
	% Bregman with BB
	opts2 = []; %modify options
    opts2.x_real = u;
    opts2.method = 0;
	tic; 
	[x2, out2] = l1l1_breg(x0, A, b, mu, opts2);
	t2 = toc;
    fprintf('Bregman with BB solved\n')
	% Bregman with Nesterov
	opts4 = []; %modify options
    opts4.x_real = u;
    opts4.method = 1;
	tic; 
	[x4, out4] = l1l1_breg(x0, A, b, mu, opts4);
	t4 = toc;
    fprintf('Bregman with Nesterov solved\n')
% 	Linearized Bregman with BB
	opts5 = []; %modify options
	tic; 
    opts5.x_real = [mu*u; zeros(m,1)];
	[x5, out5] = l1l1_lbreg_bb(x0, A, b, mu, opts5);
	t5 = toc;    
    fprintf('Linearized Bregman with BB solved\n')
	% ADMM with dual
	opts6 = []; %modify options
	tic; 
    opts6.x_real = [mu*u; zeros(m,1)];
	[x6, out6] = l1l1_dual_ADMM(x0, A, b, mu, opts6);
	t6 = toc;
    fprintf('ADMM with dual solved\n')


    
    fprintf('cvx-call-mosek: res: %3.2e, nrm1: %3.2e, cpu: %5.2f\n', resfun(x1), nrm1fun(x1), t1);
    fprintf('Bregman-BB: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x2), nrm1fun(x2), t2, errfun(x1, x2));
    fprintf('Bregman-Nesterov: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x4), nrm1fun(x4), t4, errfun(x1, x4));
    fprintf('LBreg-BB: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x5), nrm1fun(x5), t5, errfun(x1, x5));
    fprintf('ADMM-dual: res: %3.2e, nrm1: %3.2e, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', ...
        resfun(x6), nrm1fun(x6), t6, errfun(x1, x6));

    fig = figure;
    hold on
    subplot(6,1,1)
    plot(u)
    legend('solution')
    subplot(6,1,2)
    plot(x1)
    legend('cvx mosek')
    subplot(6,1,3)
    plot(x2)
    legend('Bregman with BB')
    subplot(6,1,4)
    plot(x4)
    legend('Bregman with Nesterov')
    subplot(6,1,5)
    plot(x5)
    legend('Linearized Bregman with BB')
    subplot(6,1,6)
    plot(x6)
    legend('ADMM with dual')
    hold off
end
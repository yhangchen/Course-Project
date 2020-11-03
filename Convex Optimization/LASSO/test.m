	% min 0.5 ||Ax-b||_2^2 + mu*||x||_1

    seed = 97006855;
    ss = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(ss);
    n = 1024;
    m = 512;
    A = randn(m,n);
    u = sprandn(n,1,0.1);
    b = A*u;
    mu = 1e-3;

	x0 = rand(n,1);

	errfun = @(x1, x2) norm(x1-x2)/(1+norm(x1));
    nonzero = @(x) sum(abs(x)>1e-6)/n;

	% cvx calling mosek
	opts1 = []; 
	tic; 
	[x1, out1] = l1_cvx_mosek(x0, A, b, mu, opts1);
	t1 = toc;

	% cvx calling gurobi
	opts2 = []; 
	tic; 
	[x2, out2] = l1_cvx_gurobi(x0, A, b, mu, opts2);
	t2 = toc;

	% call mosek directly
	opts3 = []; 
	tic; 
	[x3, out3] = l1_mosek(x0, A, b, mu, opts3);
	t3 = toc;

	% call gurobi directly
	opts4 = []; 
	tic; 
	[x4, out4] = l1_gurobi(x0, A, b, mu, opts4);
	t4 = toc;

	opts5 = []; 
	tic; 
	[x5, out5] = l1_PGD_primal(x0, A, b, mu, opts5);
	t5 = toc;
 
	opts6 = []; 
	tic; 
	[x6, out6] = l1_SGD_primal(x0, A, b, mu, opts6);
	t6 = toc;
    
    % gradient method for smoothed primal problem
	opts7 = []; 
	tic; 
	[x7, out7] = l1_GD_primal(x0, A, b, mu, opts7);
	t7 = toc;

	% fast gradient method for smoothed primal problem
	opts8 = []; 
	tic; 
	[x8, out8] = l1_FGD_primal(x0, A, b, mu, opts8);
	t8 = toc;
    
    % proximal gradient
	opts9 = []; 
	tic; 
	[x9, out9] = l1_ProxGD_primal(x0, A, b, mu, opts9);
	t9 = toc;

	% fast proximal gradient
	opts10 = []; 
	tic; 
	[x10, out10] = l1_FProxGD_primal(x0, A, b, mu, opts10);
	t10 = toc;
   
	% ALM dual
	opts11 = []; 
	tic; 
	[x11, out11] = l1_ALM_dual(x0, A, b, mu, opts11);
	t11 = toc;
    
    % ADMM dual
	opts12 = []; 
	tic; 
	[x12, out12] = l1_ADMM_dual(x0, A, b, mu, opts12);
	t12 = toc;
    
    % linearized ADMM for primal
	opts13 = []; 
	tic; 
	[x13, out13] = l1_ADMM_lprimal(x0, A, b, mu, opts13);
	t13 = toc;

	fprintf('    cvx-call-mosek: cpu: %5.2f, value: %3.2e, sparsity: %3.2f\n', t1, out1.val, nonzero(x1));
	fprintf('   cvx-call-gurobi: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t2, errfun(x1, x2), (out2.val-out1.val)/out1.val, nonzero(x2));
	fprintf('        call-mosek: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t3, errfun(x1, x3), (out3.val-out1.val)/out1.val, nonzero(x3));
	fprintf('       call-gurobi: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t4, errfun(x1, x4), (out4.val-out1.val)/out1.val, nonzero(x4));
	fprintf('          projgrad: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t5, errfun(x1, x5), (out5.val-out1.val)/out1.val, nonzero(x5));
	fprintf('           subgrad: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t6, errfun(x1, x6), (out6.val-out1.val)/out1.val, nonzero(x6));
    fprintf('     grad smoothed: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t7, errfun(x1, x7), (out7.val-out1.val)/out1.val, nonzero(x7));
	fprintf('fast grad smoothed: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t8, errfun(x1, x8), (out8.val-out1.val)/out1.val, nonzero(x8));
	fprintf('          proxgrad: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t9, errfun(x1, x9), (out9.val-out1.val)/out1.val, nonzero(x9));
 	fprintf('     fast proxgrad: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t10, errfun(x1, x10), (out10.val-out1.val)/out1.val, nonzero(x10));
 	fprintf('          alm dual: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t11, errfun(x1, x11), (out11.val-out1.val)/out1.val, nonzero(x11));
 	fprintf('         admm dual: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t12, errfun(x1, x12), (out12.val-out1.val)/out1.val, nonzero(x12));
  	fprintf('admm linear primal: cpu: %5.2f, err-to-cvx-mosek: %3.2e, objval to cvx mosek %3.2e, sparsity: %3.2f\n', t13, errfun(x1, x13), (out13.val-out1.val)/out1.val, nonzero(x13));
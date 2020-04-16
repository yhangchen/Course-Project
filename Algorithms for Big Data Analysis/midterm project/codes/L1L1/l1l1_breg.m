function [x,out]=l1l1_breg(x0, A, b, mu, opts)
    % implement Bregman algorithm from the paper:
    % Bregman Iterative Algorithms for \ell_1-Minimization with Applications to Compressed Sensing    
    % where \beta means \mu in the paper.
    % We use proxgrad with BB

    if ~isfield(opts,'maxIter');	opts.maxIter  = 100;  end
    if ~isfield(opts,'tol');       opts.tol       = 1e-12;   end
    if ~isfield(opts,'beta') 
        opts.beta       = 0.01;
    end
    if ~isfield(opts,'display_err'); opts.display_err       = 1;  end
    if ~isfield(opts,'method'); opts.method = 0;end
    % method 
    %       0 ista with BB
    %       1 nesterov

    maxIter = opts.maxIter;
    tol = opts.tol;

    [m0, n0] = size(A);
    xk = [mu*x0; b-A*x0];
    A = [A mu*eye(m0)]; b = mu*b;
    x_real = [mu*opts.x_real;zeros(m0,1)];
    nx_real = norm(x_real);
    display_err = opts.display_err;
    beta = opts.beta;
    
    k = 0;
    bk = b;
    err = zeros(maxIter,1);
    
    opts_in = [];
    method = opts.method;
    if method == 0
        lasso = @l1_proxgrad;
    elseif method == 1
        lasso = @l1_fprox_nesterov;
    else
        error('incorrect method');
    end
    % avoid repeat computation of ATA, set initial stepsize 1/norm(ATA)
    ATA = A'*A;
    alpha0 = 1/norm(ATA);
    opts_in.alpha = alpha0; opts_in.ATA = ATA;
    while k < maxIter
        k = k + 1;
        xk = lasso(xk, A, bk, beta, opts_in); % lasso subproblem
        if display_err
            err(k) = norm(xk-x_real)/nx_real;
            fprintf('relative error: '+string(err(k))+'\n')
        end
        resk = b - A*xk;
        bk = bk + resk;
        if norm(resk) < tol
            break
        end
    end
    if k == maxIter
        fprintf('Maximum number of Bregman iteration reached.\n');
    end
    x = xk(1:n0)/mu;
    err = err(1:k);
    out = [];
    out.iter = k;
    out.val = norm(xk, 1);
    out.err = err;
end

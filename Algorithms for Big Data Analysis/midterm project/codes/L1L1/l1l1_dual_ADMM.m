function [x, out] = l1l1_dual_ADMM(x0, A, b, mu, opts)
    % implement algorithm in section 2.3 "Applying ADM to dual problems" from the paper
    % ALTERNATING DIRECTION ALGORITHMS FOR \ell_1-PROBLEMS IN COMPRESSIVE SENSING
    % \beta, \gamma's definition and default value is the same as in the paper.
    if ~isfield(opts,'maxIter'); opts.maxIter = 1e4; end
    if ~isfield(opts,'beta'); opts.beta = 5e-3; end
    if ~isfield(opts,'gamma'); opts.gamma = (1+sqrt(5))/2; end
    if ~isfield(opts,'display_err'); opts.display_err  = 0;  end
    if ~isfield(opts,'tol'); opts.tol = 1e-10; end
    x_real = opts.x_real; nxreal = norm(x_real);
    % change l1l1 to BP
    [m0, n0] = size(A);
    A = [A mu*eye(m0)]/sqrt(1+mu^2); b = mu*b/sqrt(1+mu^2);
    
    isorth =  is_orth(A); % whether A's rows are orthogonal
    [m, n] = size(A); x = A'*b;
    z = zeros(n, 1); y = zeros(m, 1); Aty = zeros(n,1);
    normb = norm(b); maxIter = opts.maxIter; beta = opts.beta; gamma = opts.gamma; tol = opts.tol;
    for iter = 1:maxIter
        if isorth
            y = A*(z - x/beta) + b/beta;
            Aty = A'*y;
        else
            % instead of solving (AA')^{-1}£¬use steepest descent.
            g = A*(Aty - z)*beta + A*x - b;
            Atg = A'*g;
            step = norm(g)^2/(beta*norm(Atg)^2+eps);
            y = y - step*g;
            Aty = Aty - step*Atg;
        end
        z = Aty + x/beta;
        z = z ./ max(1,abs(z));
        xp = x;
        resz = z - Aty;
        x = x - gamma*beta*resz; 
        if opts.display_err
            err(iter) = norm(x-x_real)/nxreal;
            fprintf('relative error: '+string(err(iter))+'\n')
        end
        % stopping
        rel_xcge = norm(x-xp)/max(norm(x),eps); % relative change in x
        if rel_xcge > tol
            continue
        elseif rel_xcge < tol/2
            break
        end
        rel_resz = norm(resz)/max(norm(z),eps); 
        pvalue = norm(x,1);
        dvalue = b'*y;
        rel_dgap = abs(pvalue-dvalue)/max(abs(pvalue),eps); 
        rel_pgap = norm(A*x-b)/normb;
        if max([rel_resz,rel_dgap,rel_pgap]) < tol % criterion (2.32) in the reference.
            break
        end
    end
    if iter == maxIter
        fprintf('Maximum number of outer iteration reached.\n');
    end
    out.iter = iter;
    out.val = norm(x,1);
    x = x(1:n0)/mu;
end

function bool = is_orth(A)
    % check whether the rows of A are orthonormal
    bool = 0;
    s1 = randn(size(A,1),1);
    s2 = A'*s1;
    s2 = A*s2;
    err = norm(s1-s2)/norm(s1);
    if err < eps
        bool = 1;
    end
end
function x = l1_fprox_nesterov(x0, ATA, ATb, mu, opts)
% modified from L1L1 part, replace A,b with ATA, ATb
    try
        chol(ATA);
    catch
        [V, D] = eig(ATA);
        d = diag(D);
        ATA = (V * diag(max(d,0))) * V';
        ATb = V * (diag(d>0) * (V' * ATb));
    end
	if ~isfield(opts,'alpha');          opts.alpha  = 1/norm(ATA);   end
	if ~isfield(opts,'mu_ratio');   opts.mu_ratio   = 0.5;    end
	if ~isfield(opts,'maxIter');	opts.maxIter    = 1e4;	  end
	if ~isfield(opts,'tol1');       opts.tol1       = 1e-8;   end
	if ~isfield(opts,'tol2');       opts.tol2       = 1e-8;  end
   
	alpha = opts.alpha;				
	maxIter = opts.maxIter;
	mu_ratio = opts.mu_ratio;
	tol1 = opts.tol1;
	tol2 = opts.tol2;

	x = x0;
	mui = max(mu, mu_ratio*max(abs(ATb)));
	%main loop
	iter = 0;
    t = 1; xp = x;
    while mui>mu+eps && iter < maxIter
        while iter < maxIter
            y = x + (t - 1)/(t+2)*(x - xp);
            xp = x; 
            grad = ATA*y - ATb;
            y0 = y - alpha*grad;
            x = sign(y0).*max(abs(y0)-alpha*mui,0);
            iter = iter + 1;
            t = t + 1;
            if sum(abs(x-xp)) < tol1
                break
            end
        end
        mui = max(mu, mu_ratio*(min(mui, max(abs(ATA*x - ATb)))));
        t = 1; xp = x;
        y = x + (t - 1)/(t+2)*(x - xp);
        grad = ATA*y - ATb;
        y0 = y - alpha*grad;
        x = sign(y0).*max(abs(y0)-alpha*mui,0);
        iter = iter + 1;
        t = t + 1;
    end
    t = 1; xp = x;
    while iter < maxIter
        y = x + (t - 1)/(t+2)*(x - xp);
        xp = x; 
        grad = ATA*y - ATb;
        y0 = y - alpha*grad;
        x = sign(y0).*max(abs(y0)-alpha*mui,0);
        iter = iter + 1;
        t = t + 1;
        if sum(abs(x-xp)) < tol2
                break
        end
    end
end
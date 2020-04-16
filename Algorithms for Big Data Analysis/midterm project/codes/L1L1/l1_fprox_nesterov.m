function [x,out] = l1_fprox_nesterov(x0, A, b, mu, opts)
	if ~isfield(opts,'alpha');          opts.alpha  = 5e-4;   end
	if ~isfield(opts,'mu_ratio');   opts.mu_ratio   = 0.5;    end
	if ~isfield(opts,'maxIter');	opts.maxIter    = 1e4;	  end
	if ~isfield(opts,'tol1');       opts.tol1       = 1e-6;   end
	if ~isfield(opts,'tol2');       opts.tol2       = 1e-6;  end
    if ~isfield(opts,'ATA');    	    opts.ATA       	= A'*A;   	  end
    if ~isfield(opts,'ATb');    	    opts.ATb       	= A'*b;   	  end

	alpha = opts.alpha;				
	maxIter = opts.maxIter;
	mu_ratio = opts.mu_ratio;
	tol1 = opts.tol1;
	tol2 = opts.tol2;
    ATA = opts.ATA;
    ATb = opts.ATb;
	x = x0;
	mui = max(mu, mu_ratio*max(abs(ATb)));
	f = 0.5*norm(A*x-b)^2+mui*norm(x,1);
	ratio = 1;
	%main loop
	iter = 0;
    t = 1; xp = x;
    while mui>mu+eps && iter < maxIter
        while ratio>=tol1 && iter < maxIter
            y = x + (t - 1)/(t+2)*(x - xp);
            xp = x; 
            grad = ATA*y - ATb;
            y0 = y - alpha*grad;
            x = sign(y0).*max(abs(y0)-alpha*mui,0);
            fp = f; 
            f = 0.5*norm(A*x-b)^2+mui*norm(x,1);
            ratio = abs((f-fp)/fp);
            iter = iter + 1;
            t = t + 1;
        end
        mui = max(mu, mu_ratio*(min(mui, max(abs(ATA*x - ATb)))));
        t = 1; xp = x;
        y = x + (t - 1)/(t+2)*(x - xp);
        grad = ATA*y - ATb;
        y0 = y - alpha*grad;
        x = sign(y0).*max(abs(y0)-alpha*mui,0);
        fp = f; 
        f = 0.5*norm(A*x-b)^2+mui*norm(x,1);
        ratio = abs((f-fp)/fp);
        iter = iter + 1;
        t = t + 1;
    end
    t = 1; xp = x;
    while ratio>=tol2 && iter < maxIter
        y = x + (t - 1)/(t+2)*(x - xp);
        xp = x; 
        grad = ATA*y - ATb;
        y0 = y - alpha*grad;
        x = sign(y0).*max(abs(y0)-alpha*mui,0);
        fp = f; 
        f = 0.5*norm(A*x-b)^2+mui*norm(x,1);
        ratio = abs((f-fp)/fp);
        iter = iter + 1;
        t = t + 1;
    end
    out.val = 0.5*norm(A*x-b)^2+mui*norm(x,1);
end
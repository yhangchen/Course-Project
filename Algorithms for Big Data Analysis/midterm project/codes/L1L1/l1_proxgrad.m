function [x,out] = l1_proxgrad(x0, A, b, mu, opts)
    %This program implements the proximal gradient 
    %method with continuation method.

    if ~isfield(opts,'alpha');          opts.alpha  = 5e-4;   end
    if ~isfield(opts,'mu_ratio'); opts.mu_ratio = 0.5;        end
    if ~isfield(opts,'maxIter');	opts.maxIter    = 1e4;	  end
    if ~isfield(opts,'tol1');       opts.tol1       = 1e-6;   end
    if ~isfield(opts,'tol2');       opts.tol2       = 1e-6;  end % setting this to 1e-5 will be faster in rng(1234)
    if ~isfield(opts,'BB');    	    opts.BB       	= 1;   	  end
    if ~isfield(opts,'ATA');    	    opts.ATA       	= A'*A;   	  end
    if ~isfield(opts,'ATb');    	    opts.ATb       	= A'*b;   	  end
    
    alpha = opts.alpha;				
    maxIter = opts.maxIter;
    mu_ratio = opts.mu_ratio;
    tol1 = opts.tol1;
    tol2 = opts.tol2;

    ATA = opts.ATA;
    ATb = opts.ATb;

    x = x0; xp = x;
    mui = max(mu, mu_ratio* max(abs(ATb(:))));
    grad = ATA*x - ATb;
    x0 = x - alpha*grad;
    x = sign(x0).*max(abs(x0)-alpha*mui,0);
    gradp = grad;
    grad = ATA*x - ATb;
    f = 0.5*norm(A*x-b)^2+ mui*norm(x,1);
    ratio = 1;
    %main loop
    iter=0;
    while mui>mu+eps && iter < maxIter
        while ratio>=tol1 && iter < maxIter
            if opts.BB
                dx = x - xp;
                dgrad = grad -gradp;
                alpha = dx'*dx/(dgrad'*dx);
            end
            fp = f; xp = x; gradp = grad;
            x0 = x - alpha*grad;
            x = sign(x0).*max(abs(x0)-alpha*mui,0);
            grad = ATA*x - ATb;
            f = 0.5*norm(A*x-b)^2 + mui*norm(x,1);
            ratio = abs((f-fp)/fp);
            iter = iter+1;
        end
        mui = max(mu, mu_ratio*(min(mui, max(abs(grad(:))))));
        alpha = opts.alpha;
        x0 = x - alpha*grad;
        x = sign(x0).*max(abs(x0)-alpha*mui,0);
        grad = ATA*x - ATb;
        f = 0.5*norm(A*x-b)^2 + mui*norm(x,1);
        ratio = abs((f-fp)/fp);
        iter = iter+1;
    end

    while ratio>=tol2 && iter < maxIter
        if opts.BB
            dx = x - xp;
            dgrad = grad -gradp;
            alpha = dx'*dx/(dgrad'*dx);
        end
        fp = f; xp = x; gradp = grad;
        x0 = x - alpha*grad;
        x = sign(x0).*max(abs(x0)-alpha*mui,0);
        grad = ATA*x - ATb;
        f = 0.5*norm(A*x-b)^2 + mui*norm(x,1);
        ratio = abs(f-fp)/max(1,abs(fp));
        iter = iter+1;
    end
    out = [];
    out.val = 0.5*norm(A*x-b)^2+mu*norm(x,1);
end
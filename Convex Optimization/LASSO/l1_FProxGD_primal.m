function [x,out] = l1_FProxGD_primal(x0, A, b, mu, opts)
	if ~isfield(opts,'alpha');      opts.alpha      = 4e-4;   end
	if ~isfield(opts,'gamma');      opts.gamma      = 0.5;    end
	if ~isfield(opts,'maxIter');	opts.maxIter    = 1e4;	  end
	if ~isfield(opts,'M1');         opts.M1         = 60;     end
	if ~isfield(opts,'M2');         opts.M2         = 200;     end

	alpha = opts.alpha;				
	maxIter = opts.maxIter;
	gamma = opts.gamma;

	ATA = A'*A;
	ATb = A'*b;
	bTb = b'*b;

	x = x0;
	mui = max(mu, gamma * max(abs(ATb(:))));
    grad = ATA*x - ATb;
	f = 0.5*((grad-ATb)'*x + bTb) + mui*norm(x,1);
	%main loop
	k = 1;iter = 1;inner_iter = 0;
	while mui>mu
		while inner_iter < opts.M1
			xp = x;
			y = x + (k-1)/(k+2)*(x-xp);
			grad = ATA*y - ATb;
			x = soft(y - alpha*grad, alpha*mui);
			k = k+1;
			iter = iter+1;inner_iter = inner_iter + 1;
			if (iter>maxIter);break;end
		end
		mui = max(mu, gamma*(min(mui, max(abs(ATA*x - ATb)))));
		fp = f; xp = x;
		k = 1;
        y = x + (k-1)/(k+2)*(x-xp);
        grad = ATA*y - ATb;
        x = soft(y - alpha*grad, alpha*mui);
		k = k+1; iter = iter+1;inner_iter = 0;
		if (iter>maxIter);break;end
    end

    while inner_iter < opts.M2
        xp = x;
        y = x + (k-1)/(k+2)*(x-xp);
        grad = ATA*y - ATb;
        x = soft(y - alpha*grad, alpha*mui);
        iter = iter+1;
        k = k+1;inner_iter = inner_iter + 1;
        if (iter>maxIter);break;end
    end

    out.iter = iter;
	out.val = 0.5*((grad-ATb)'*x + bTb) + mui*norm(x,1);
end

function s = soft(x,mu)
    s = sign(x).*max(abs(x)-mu,0);
end
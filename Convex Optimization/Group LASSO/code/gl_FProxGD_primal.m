function [x,iter,out] = gl_FProxGD_primal(x0, A, b, mu, opts)
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
	%main loop
	k = 1;iter = 1;inner_iter = 0;
	while mui>mu
		while inner_iter < opts.M1
			xp = x;
			y = x + (k-1)/(k+2)*(x-xp);
			grad = ATA*y - ATb;
			x = prox(y - alpha*grad, alpha*mui);
			k = k+1;
			iter = iter+1;inner_iter = inner_iter + 1;
			if (iter>maxIter);break;end
		end
		mui = max(mu, gamma*(min(mui, max(abs(grad(:))))));
		xp = x;
		k = 1;
        y = x + (k-1)/(k+2)*(x-xp);
        grad = ATA*y - ATb;
        x = prox(y - alpha*grad, alpha*mui);
		k = k+1; iter = iter+1;inner_iter = 0;
		if (iter>maxIter);break;end
    end

    while inner_iter < opts.M2
        xp = x;
        y = x + (k-1)/(k+2)*(x-xp);
        grad = ATA*y - ATb;
        x = prox(y - alpha*grad, alpha*mui);
        iter = iter+1; k = k+1;inner_iter = inner_iter + 1;
        if (iter>maxIter);break;end
    end
	out.val = 0.5*norm(A*x-b,'fro')^2+mui*sum(norms(x,2,2));
end

function s = prox(x,mu)
    s = zeros(size(x));
    for i = 1:size(x,1)
        xi = x(i,:);
        norm_xi = norm(xi);
        if norm_xi > mu
            norm_si = norm(xi) - mu;
            s(i,:) = x(i,:) * norm_si/norm_xi;
        end
    end
end
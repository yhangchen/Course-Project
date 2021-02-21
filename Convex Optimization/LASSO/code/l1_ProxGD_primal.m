function [x,out] = l1_ProxGD_primal(x0, A, b, mu, opts)
	if ~isfield(opts,'alpha');      opts.alpha      = 4e-4;   end
	if ~isfield(opts,'gamma');      opts.gamma      = 0.5;    end
	if ~isfield(opts,'maxIter');	opts.maxIter    = 1e4;	  end
	if ~isfield(opts,'M1');         opts.M1         = 20;     end
	if ~isfield(opts,'M2');         opts.M2         = 40;     end
	if ~isfield(opts,'BB');    	    opts.BB       	= 1;   	  end

	alpha = opts.alpha;				
	maxIter = opts.maxIter;
	gamma = opts.gamma;

	ATA = A'*A;
	ATb = A'*b;
	bTb = b'*b;

	x = x0;
	mui = max(mu, gamma * max(abs(ATb(:))));
	grad = ATA*x - ATb;
	gradp = grad;
	f = 0.5*((grad-ATb)'*x + bTb) + mui*norm(x,1);
	%main loop
	iter = 1;
    inner_iter = 0;
    while mui>mu
        while inner_iter < opts.M1
            xp = x;
            x = soft(x - alpha*grad, alpha*mui);
            grad = ATA*x - ATb;
            if opts.BB
                dx = x - xp;
                dgrad = grad -gradp;
                alpha = dx'*dx/(dgrad'*dx);
            end
            gradp = grad;
            iter = iter + 1; inner_iter = inner_iter + 1;
            if (iter > maxIter);break;end
        end
        mui = max(mu, gamma*(min(mui, max(abs(grad(:))))));
        alpha = opts.alpha;
        grad = ATA*x - ATb;
        f = 0.5*((grad-ATb)'*x + bTb) + mui*norm(x,1);
        iter = iter+1; inner_iter = 0;
        if (iter > maxIter);break;end
    end

    while inner_iter < opts.M2
        xp = x;
        x = soft(x - alpha*grad, alpha*mui);
        grad = ATA*x - ATb;
        if opts.BB
            dx = x - xp;
            dgrad = grad -gradp;
            alpha = dx'*dx/(dgrad'*dx);
        end
        gradp = grad;
        iter = iter+1; inner_iter = inner_iter + 1;
        if (iter > maxIter);break;end
    end
    out.iter = iter;
	out.val = 0.5*norm(A*x-b)^2+mui*norm(x,1);
end

function s = soft(x,mu)
    s = sign(x).*max(abs(x)-mu,0);
end
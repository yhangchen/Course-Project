function [x,iter, out] = gl_ProxGD_primal(x0, A, b, mu, opts)
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
	%main loop
	iter = 1;
    inner_iter = 0;
    while mui>mu
        while inner_iter < opts.M1
            xp = x;
            x = prox(x - alpha*grad, alpha*mui);
            grad = ATA*x - ATb;
            if opts.BB
                dx = x - xp;
                dgrad = grad -gradp;
                alpha = vec(dx)'*vec(dx)/(vec(dgrad)'*vec(dx));
            end
            gradp = grad;
            iter = iter + 1; inner_iter = inner_iter + 1;
            if (iter > maxIter);break;end
        end
        mui = max(mu, gamma*(min(mui, max(abs(vec(ATA*x-ATb))))));
        alpha = opts.alpha;
        grad = ATA*x - ATb;
        iter = iter+1; inner_iter = 0;
        if (iter > maxIter);break;end
    end

    while inner_iter < opts.M2
        xp = x;
        x = prox(x - alpha*grad, alpha*mui);
        grad = ATA*x - ATb;
        if opts.BB
            dx = x - xp;
            dgrad = grad -gradp;
            alpha = vec(dx)'*vec(dx)/(vec(dgrad)'*vec(dx));
        end
        gradp = grad;
        iter = iter+1; inner_iter = inner_iter + 1;
        if (iter > maxIter);break;end
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
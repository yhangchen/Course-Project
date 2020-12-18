function [x,iter, out] = gl_GD_primal(x0, A, b, mu, opts)
    if ~isfield(opts,'alpha');      opts.alpha      = 4e-4;   end
	if ~isfield(opts,'gamma');      opts.gamma      = 0.1;    end
	if ~isfield(opts,'lambda');    	opts.lambda    	= 1e-6;   end
	if ~isfield(opts,'lambdaj');    opts.lambdaj    = 1e-1;   end
	if ~isfield(opts,'beta');       opts.beta       = 0.1;    end
	if ~isfield(opts,'M_1');	    opts.M_1        = 200;	  end
	if ~isfield(opts,'M_2');	    opts.M_2        = 800;	  end
	if ~isfield(opts,'BB');    	    opts.BB         = 1;   	  end
    %copy paramter
	alpha = opts.alpha;				
	gamma = opts.gamma;
	lambdaj = opts.lambdaj;
	lambda = opts.lambda;
	beta = opts.beta;
	M_1 = opts.M_1;
	M_2 = opts.M_2;
	ATA = A'*A;
	ATb = A'*b;
	x = x0;
	mui = max(mu, gamma* max(abs(ATb(:))));
	grad = ATA*x-ATb + mui*groupgrad(x, lambdaj);
	xp = x; gradp = grad;
	x = x - alpha*grad;
	grad = ATA*x-ATb + mui*groupgrad(x, lambdaj);
    k=0;
    while (lambdaj > lambda) || (mui>mu)
        iter = 0;
        while iter<M_1-1
            if opts.BB
                dx = x - xp;
                dgrad = grad -gradp;
                alpha = vec(dx)'*vec(dx)/(vec(dgrad)'*vec(dx)+eps);
            end
            gradp = grad; xp = x;
            x = x - alpha*grad;
            grad = ATA*x-ATb + mui*groupgrad(x, lambdaj);
            k = k+1;
            iter = iter+1;
        end
        mui = max(mu, gamma*(min(mui, max(abs(vec(ATA*x-ATb))))));
        lambdaj = max(lambdaj*beta, lambda);
        alpha = min(opts.alpha,lambdaj);
        grad = ATA*x-ATb + mui*groupgrad(x, lambdaj);
        xp = x; gradp = grad;
        x = x - alpha*grad;
        grad = ATA*x-ATb + mui*groupgrad(x, lambdaj);
        k = k+1;
    end
	iter = 0;
    while iter<M_2
        if opts.BB
            dx = x - xp;
            dgrad = grad -gradp;
            alpha = vec(dx)'*vec(dx)/(vec(dgrad)'*vec(dx)+eps);
        end
        gradp = grad; xp = x;
        x = x - alpha*grad;
        grad = ATA*x-ATb + mui*groupgrad(x, lambdaj);
        k = k+1;
        iter = iter+1;
    end
    iter = k;
	out.val = 0.5*norm(A*x-b,'fro')^2+mu*sum(norms(x,2,2));
end


function x = groupgrad(x, mu)
    for i = 1:size(x,1)
        if norm(x(i,:)) > mu
            x(i,:) = x(i,:)/norm(x(i,:));
        else
            x(i,:) = x(i,:)/mu;
        end
    end
end
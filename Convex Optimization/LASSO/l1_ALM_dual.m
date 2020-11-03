function [x,out] = l1_ALM_dual(x0, A, b, mu, opts)
	if ~isfield(opts,'gamma');      opts.gamma      = 0.1;    end
	if ~isfield(opts,'maxIter');	opts.maxIter    = 1e4;	  end
    if ~isfield(opts,'M1');         opts.M1         = 8;     end
    if ~isfield(opts,'M2');         opts.M2         = 8;     end
    if ~isfield(opts,'M3');         opts.M3         = 4;      end
    if ~isfield(opts,'a');          opts.a          = 5e-2;   end
    a = opts.a;
    gamma = opts.gamma;
    maxIter = opts.maxIter;
    L = inf;
    [m, n] = size(A);
%     ATA = A'*A;
    ATb = A'*b;
    iter = 0;
	mui = max(mu, gamma * max(abs(ATb(:))));
    z = zeros(m,1); lambda = zeros(n,1);
    H = eye(m) + a*(A*A');
    [LH,DH] = ldl(H);
    outer_iter = 0;
    while mui>mu
		while outer_iter < opts.M1
            inner_iter = 0;
            while inner_iter < opts.M3
                [d] = newton(z, lambda);
                z = z - LH'\(DH\(LH\d));
                inner_iter = inner_iter + 1;
            end
			v = lambda/a+A'*z;
			w = v - soft(v,mui);
			lambda = lambda +a *(v-lambda/a-w);
			iter = iter + 1;
            outer_iter = outer_iter + 1;
			if (iter>maxIter);break;end
		end
		mui = max(mu, gamma*mui);
		iter = iter+1;
        outer_iter = 0;
		if (iter>maxIter);break;end
    end
    
    while outer_iter < opts.M2
        inner_iter = 0;
        while inner_iter < opts.M3
            [d] = newton(z, lambda);
            z = z - LH'\(DH\(LH\d));
            inner_iter = inner_iter + 1;
        end
        v = lambda/a+A'*z;
        w = v - soft(v,mui);
        lambda = lambda + a*(v-lambda/a-w);
        iter = iter + 1;
        outer_iter = outer_iter + 1;
        if (iter>maxIter);break;end
    end
    
    x = -lambda;
    out.iter = iter;
	out.val = 0.5*norm(A*x-b)^2+mui*norm(x,1);

    function [d] = newton(z, lambda)
        v = lambda/a+A'*z;
        d = z+b+a*A*soft(v,mui);
    end

end

function s = soft(x,mu)
    s = sign(x).*max(abs(x)-mu,0);
end

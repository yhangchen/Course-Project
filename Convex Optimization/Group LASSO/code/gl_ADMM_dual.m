function [x,iter,out] = gl_ADMM_dual(x0, A, b, mu, opts)
	if ~isfield(opts,'gamma');      opts.gamma      = 0.2;    end
    if ~isfield(opts,'M1');         opts.M1         = 10;     end
    if ~isfield(opts,'M2');         opts.M2         = 15;      end
    if ~isfield(opts,'a');          opts.a          = 1e-2;   end
    a = opts.a;
    gamma = opts.gamma;
    [m, n] = size(A);
    ATb = A'*b;
    iter = 0;
	mui = max(mu, gamma * max(abs(ATb(:))));
    z = zeros(m,1); lambda = zeros(n,1);
    H = eye(m) + a*(A*A');    
    [LH,DH] = ldl(H);
    inner_iter = 0;
    while mui>mu
		while inner_iter < opts.M1
			v = lambda/a + A'*z;
			w = v - prox(v,mui);
			z = LH'\(DH\(LH\(-b-A*(lambda-a*w))));
            lambda = lambda + a*(A'*z-w);
			iter = iter + 1;
            inner_iter = inner_iter + 1;
		end
		mui = max(mu, gamma*mui);
        inner_iter = 0;
    end
    while inner_iter < opts.M2
        v = lambda/a + A'*z;
        w = v - prox(v,mui);
        z = LH'\(DH\(LH\(-b-A*(lambda-a*w))));
        lambda = lambda + a*(A'*z-w);
        iter = iter + 1;
        inner_iter = inner_iter + 1;
    end
    x = -lambda;
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

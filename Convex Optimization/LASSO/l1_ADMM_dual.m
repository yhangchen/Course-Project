function [x,out] = l1_ADMM_dual(x0, A, b, mu, opts)
	if ~isfield(opts,'gamma');      opts.gamma      = 0.2;    end
    if ~isfield(opts,'M1');         opts.M1         = 6;     end
    if ~isfield(opts,'M2');         opts.M2         = 20;      end
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
			w = v - soft(v,mui);
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
        w = v - soft(v,mui);
        z = LH'\(DH\(LH\(-b-A*(lambda-a*w))));
        lambda = lambda + a*(A'*z-w);
        iter = iter + 1;
        inner_iter = inner_iter + 1;
    end
    x = -lambda;
    out.iter = iter;
	out.val = 0.5*norm(A*x-b)^2+mui*norm(x,1);
end

function s = soft(x,mu)
    s = sign(x).*max(abs(x)-mu,0);
end


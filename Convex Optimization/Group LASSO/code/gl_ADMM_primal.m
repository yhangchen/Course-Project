function [x,iter,out] = gl_ADMM_primal(x0, A, b, mu, opts)
	if ~isfield(opts,'gamma');      opts.gamma      = 0.2;    end
    if ~isfield(opts,'M1');         opts.M1         = 10;     end
    if ~isfield(opts,'M2');         opts.M2         = 15;      end
    if ~isfield(opts,'a');          opts.a          = 1e2;   end
    a = opts.a;
    gamma = opts.gamma;
    [m, n] = size(A);
    ATA = A'*A;
    ATb = A'*b;
    iter = 0;
	mui = max(mu, gamma * max(abs(ATb(:))));
    H = a * eye(n) + (A'*A);    
    [LH,DH] = ldl(H);
	x = x0;	y = x0;	z = zeros(n,1); inner_iter = 0;
    while mui>mu
        while inner_iter < opts.M1
            x = LH'\(DH\(LH\(ATb+a*y-z)));
            y = prox(x+z/a,mui/a);
            z = z+a*(x-y);
            inner_iter = inner_iter + 1;
            iter = iter + 1;
        end
        mui = max(mu, gamma*(min(mui, max(vec(ATA*x-ATb))))); inner_iter = 0;
    end

    while inner_iter < opts.M2
        x = LH'\(DH\(LH\(ATb+a*y-z)));
        y = prox(x+z/a,mui/a);
        z = z+a*(x-y);
        inner_iter = inner_iter + 1;
        iter = iter + 1;
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

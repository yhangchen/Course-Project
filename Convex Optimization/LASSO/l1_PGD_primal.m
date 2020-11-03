function [x,out] = l1_PGD_primal(x0, A, b, mu, opts)
	if ~isfield(opts,'alpha');      opts.alpha       = 4e-4;   end
	if ~isfield(opts,'maxIter');    opts.maxIter     = 200;    end
	if ~isfield(opts,'numConti');   opts.numConti    = 6;      end
	if ~isfield(opts,'gamma');      opts.gamma       = 10;     end
    alpha = opts.alpha;
	maxIter = opts.maxIter;
	numConti = opts.numConti;
	gamma = opts.gamma;
    [~, n] = size(A);
	ATA = A'*A;
	ATb = A'*b;
	mui = mu*gamma^(numConti-1);
	xp = x0.*(x0>0);
	xn = xp - x0;
	%main loop
	for i=1:numConti
		for j=1:maxIter
			new_x_p = xp - alpha*(ATA*(xp-xn)+mui*ones(n,1)-ATb);
			new_x_n = xn - alpha*(ATA*(xn-xp)+mui*ones(n,1)+ATb);
			xp = new_x_p.*(new_x_p>0);
			xn = new_x_n.*(new_x_n>0);
		end
		if i<numConti
			mui = mui/gamma;
		end
    end
	x = xp-xn;
	out = [];
	out.val = 0.5*norm(A*x-b)^2+mu*norm(x,1);
    out.alpha = alpha;
	out.mu = mui;
end
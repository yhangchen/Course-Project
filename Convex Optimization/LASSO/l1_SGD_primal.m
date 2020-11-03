function [x,out] = l1_SGD_primal(x0, A, b, mu, opts)
	if ~isfield(opts,'alpha');      opts.alpha      = 4e-4; end
	if ~isfield(opts,'maxIter');    opts.maxIter    = 400;    end
	if ~isfield(opts,'numConti');   opts.numConti   = 6;      end
	if ~isfield(opts,'gamma');      opts.gamma      = 10;     end
	alpha = opts.alpha;		
	maxIter = opts.maxIter;
	numConti = opts.numConti;
	gamma = opts.gamma;
	ATA=A'*A;
	ATb=A'*b;
	mui = mu*gamma^(numConti-1);
	x = x0;
    for i=1:numConti
        for j=1:maxIter
            x=x-alpha*(ATA*x-ATb+mui*sign(x));
        end
        if i<numConti
            mui = mui/gamma;
        end
    end
	out.val = 0.5*norm(A*x-b)^2+mui*norm(x,1);
	out.alpha = alpha;
	out.mu = mui;
end
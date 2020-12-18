function [x,iter,out] = gl_SGD_primal(x0, A, b, mu, opts)
if ~isfield(opts,'alpha');      opts.alpha       = 4e-4;  end
if ~isfield(opts,'maxIter');    opts.maxIter     = 200;  end
if ~isfield(opts,'numConti');   opts.numhomotopy = 10;    end
if ~isfield(opts,'gamma');      opts.gamma       = 4;     end
alpha = opts.alpha;
maxIter = opts.maxIter;
numhomotopy = opts.numhomotopy;
gamma = opts.gamma;
ATA=A'*A;
ATb=A'*b;
mui = mu*gamma^(numhomotopy-1);
x = x0;
l = size(x,2);
iter = 0;
for i=1:numhomotopy
    if i == numhomotopy
        maxIter = 1000;
    end
    for j=1:maxIter
        x=x-alpha*(ATA*x-ATb+mui*(x.*((x.*x)*ones(l,l)).^(-1/2)));
        iter = iter + 1;
    end
    if i<numhomotopy
        mui = mui/gamma;
    end
end
out.val = 0.5*norm(A*x-b,'fro')^2+mui*sum(norms(x,2,2));
out.alpha = alpha;
out.mu = mui;
end
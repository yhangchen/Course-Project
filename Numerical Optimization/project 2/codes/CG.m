function [x, out] = CG(fun, x0, method, opts)
%{
    Non-linear conjugate gradient methods

    Parameters
    ----------
    fun: 
        objective function, with callable method f and g
    x0: 
        initial point
    method: 
        options: 'FR' for FR, 'PRP' for PRP, 'PRP+' for PRP+, 'HS' for HS,
                 'CD' for conjugate descent, 'DY' for Dai-Yuan, 'FR-PRP'
                 for FR-PRP
    opts.linesearch:
        'exact' for exact line search, 'inexact' for inexact line search
    opts.tol: 
        tolerance, used for stopping criterion, default 1e-6
    opts.maxIter
        maximum number of iterations, default 1e4
    opts.n_restart: 
        restart the iteration every n_restart steps
    opts.debug: 
        output information for every iteration if set to 1

    Returns
    -------
    x: 
        optimal point
    out.value: 
        optimal function value
    out.gradnorm: 
        norm of gradient at optimal point
    out.iter: 
        number of iterations
    out.feval:
        number of function evaluations (f and g)
    out.restart:
        number of restarting the iteration.

%}
if nargin < 4; opts = []; end
if ~isfield(opts,'linesearch'); opts.linesearch = 'inexact'; end
if ~isfield(opts,'tol'); opts.tol = 1e-6; end
if ~isfield(opts,'maxIter'); opts.maxIter = 1e4; end
if ~isfield(opts,'n_restart'); opts.n_restart = opts.maxIter; end % without periodic restart
if ~isfield(opts,'debug'); opts.debug = 1; end
x = x0;
n = length(x);
f0 = -inf; g0 = zeros(n,1);
[f1, g1] = fun(x);
d = -g1;
iter = 0; feval = 2; restart = 0; alpha = 1;
nf1 = max(abs(f1),1); ng1 = max(norm(g1),1);
xs = {}; values = []; gradnorms = [];
while norm(g1,inf) > opts.tol*(1+abs(f1))
    if abs(f1-f0)<eps
        break
    end
    if dot(d,g1) > eps || ~mod(iter,opts.n_restart)
            d = -g1;
            restart = restart + 1;
            if opts.debug
                fprintf('Restart.\n');
            end
    end
    iter = iter + 1;
    alpha0 = alpha;
    try
        [alpha, info] = bolinesearch(fun, x, d);
    catch 
        info(1) = 1;
    end
    if info(1) % stepsize is not properly selected.
        Rule.opt = [0 10 25 eps];
        [alpha, info] = bolinesearch(fun, x, d, Rule);
        if ~info(1) || abs(alpha) < eps
            alpha = alpha0;
        end
    end
    x = x + alpha*d;
    f0 = f1; g0 = g1;
    [f1, g1] = fun(x);
    feval = feval + info(3) + 2;
    if strcmp(method,'FR'); beta = dot(g1,g1)/dot(g0,g0);
    elseif strcmp(method,'PRP'); beta = dot(g1,g1-g0)/dot(g0,g0); 
    elseif strcmp(method,'PRP+'); beta = max(dot(g1,g1-g0)/dot(g0,g0),0);
    elseif strcmp(method,'CD'); beta = -dot(g1,g1)/dot(g0,d);
    elseif strcmp(method,'DY'); beta = dot(g1,g1)/dot(g1-g0,d);
    elseif strcmp(method,'HS'); beta = dot(g1,g1-g0)/dot(g1-g0,d);
    elseif strcmp(method, 'FR-PRP')
        beta1 = dot(g1,g1-g0)/dot(g0,g0); % PRP
        beta2 = dot(g1,g1)/dot(g0,g0); % FR
        beta = min(beta2,max(-beta2,beta1));
    else; error('Invalid method!\n'); 
    end
    d = beta*d - g1;
    if opts.debug
        fprintf(string(method)+' iter: '+string(iter)+', alpha: '+ string(alpha) + ...
            ', beta: '+string(beta)+ ', value: ' + string(f1)+', gradnorm: '+ string(norm(g1,inf))+...
            ', directnorm: '+ string(norm(d))+ '.\n');
        xs{iter} = x;
        values = [values;f1];
        gradnorms = [gradnorms;norm(g1,inf)];
    end
    if iter == opts.maxIter
        fprintf('Max iteration achieved, the solution might not be optimal.\n');
        break
    end
end
out = [];
out.value = f1;
out.gradnorm = norm(g1,inf);
out.iter = iter;
out.feval = feval;
out.restart = restart;
out.xs = xs;
out.values = values;
out.gradnorms = gradnorms;
end


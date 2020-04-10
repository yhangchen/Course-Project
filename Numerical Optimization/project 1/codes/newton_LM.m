function [x, out] = newton_LM(fun, x0, opts)
%{
    Modified Newton's method: LM method

    Parameters
    ----------
    fun: objective function, with callable method f, g and G
    x0: initial point
    opts.linesearch: optional
        'exact' for exact line search, 'inexact' for inexact line search
    opts.tol1: optional
        tolerance for convergence criterion
    opts.tol2: optional
        tolerance for convergence criterion

    opts.maxiter: optional
        maximum number of iterations

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
    out.feva:
        number of function evaluations
    """
%}
if ~isfield(opts,'maxiter');	opts.maxiter    = 1e3;       end
if ~isfield(opts,'tol1');       opts.tol1       = 1e-8;      end
if ~isfield(opts,'tol2');       opts.tol2       = 1e-8;      end
if ~isfield(opts,'linesearch'); opts.linesearch = 'inexact'; end
if ~isfield(opts,'dnorm');      opts.dnorm      = 1;         end


fp = -inf;
[fk, gk, Gk] = fun(x0);
xk = x0;
iter = 0;
feva = 0;
while (abs(fp-fk) > opts.tol1 && norm(gk) > opts.tol2) && iter < opts.maxiter
    vk = 0;
    while true
       try 
           R = chol(Gk + vk * eye(size(Gk)));
            break
       catch ME
            if vk == 0
                vk = norm(Gk)/2;
            else
                vk = vk * 2;
            end
       end
    end
    dk = - R'\gk;
    dk = R\dk;
    dk = opts.dnorm*dk/norm(dk);
    [alpha, info] = bolinesearch(fun, xk, dk);
    ineva = info(3);
    xk = xk + alpha * dk;
    fp = fk;
    [fk, gk, Gk] = fun(xk);
    iter = iter + 1;
    feva = feva + ineva + 3;
end
x = xk;
out = [];
out.value = fk;
out.gradnorm = norm(gk);
out.iter = iter;
out.feva = feva;
end
function [x, out] = newton_vanilla(fun, x0, opts)
%{
    Vanilla Newton's method
    Parameters
    ----------
    fun: objective function, return f, g and G
    x0: initial point
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
fp = -inf;
[fk, gk, Gk] = fun(x0);
xk = x0;
iter = 0;
feva = 3;
while (abs(fp-fk) > opts.tol1 && norm(gk) > opts.tol2) && iter < opts.maxiter
    try 
        chol(Gk);
    catch ME
       error('Hessian is not positive definite')
    end
    dk = -Gk\gk;
    xk = xk + dk;
    fp = fk;
    [fk, gk, Gk] = fun(xk);
    iter = iter + 1;
    feva = feva + 3;
end
x = xk;
out = [];
out.value = fk;
out.gradnorm = norm(gk);
out.iter = iter;
out.feva = feva;
end
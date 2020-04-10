function [x, out] = newton_mixed(fun, x0, opts)
%{
    Modified Newton's method: mixed direction

    Parameters
    ----------
    fun: objective function, with callable method f, g and G
    x0: initial point
    opts.linesearch: optional
        'exact' for exa    opts.tol2: optional
        tolerance for convergence criterion
    opts.tol3: optional
        tolerance for being orthogonal

    opts.maxiter: optional
        maximum number of iterations
ct line search, 'inexact' for inexact line search
    opts.tol1: optional
        tolerance for convergence criterion

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
if ~isfield(opts,'tol1');       opts.tol1       = 1e-12;      end
if ~isfield(opts,'tol2');       opts.tol2       = 1e-12;      end
if ~isfield(opts,'tol3');       opts.tol3       = 1e-12;      end
if ~isfield(opts,'tol4');       opts.tol4       = 1e-12;      end
if ~isfield(opts,'linesearch'); opts.linesearch = 'inexact'; end
if ~isfield(opts,'dnorm');      opts.dnorm      = 1;         end


fp = -inf;
[fk, gk, Gk] = fun(x0);
xk = x0;
iter = 0;
feva = 3;
while (abs(fp-fk) > opts.tol1 && norm(gk) > opts.tol2) && iter < opts.maxiter
    if cond(Gk) > 1e15
        cond(Gk)
        % to find out if G is singular.
        % note that det is not an appropriate choice
        % since e.g. det(0.0001*eye(100))=0
        dk = -gk;
    else
        dk = -Gk\gk;
        if abs(gk'*dk) < opts.tol3 * norm(gk) * norm(dk)
            dk = -gk;
        elseif gk'*dk > opts.tol4 * norm(gk) * norm(dk)
            dk = -dk;
        end
    end
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
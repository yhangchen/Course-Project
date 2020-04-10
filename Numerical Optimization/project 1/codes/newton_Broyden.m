function [x, out] = newton_Broyden(fun, x0, phi, opts)
%{
    Modified Newton's method: Broyden class method

    Parameters
    ----------
    fun: objective function, with callable method f, g and G
    x0: initial point
    phi: parameter of Broyden method, in [0,1]
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
[fk, gk] = fun(x0);
n = length(gk);
xk = x0;
iter = 0;
feva = 2;
Hk = eye(n);
while (abs(fp-fk) > opts.tol1 && norm(gk) > opts.tol2) && iter < opts.maxiter
    dk = -Hk * gk;
    dk = opts.dnorm*dk/norm(dk);
    try
    [alpha, info] = bolinesearch(fun, xk, dk);
    ineva = info(3);
    catch
    alpha = 1; ineva = 0;
    end
    xk = xk + alpha * dk;
    fp = fk; gp = gk;
    [fk, gk] = fun(xk);
    sk = alpha * dk; yk = gk - gp;
    Hyk = Hk*yk; syk = sk'*yk; yHyk = yk'*Hyk;
    if abs(syk) < eps 
        syk = eps;
    end
    if abs(yHyk) < eps
        yHyk = eps;
    end
    vk = sqrt(yHyk)*sk/syk-Hyk/sqrt(yHyk);
    Hk = Hk + (sk*sk')/syk - (Hyk*Hyk')/yHyk + phi * (vk*vk');
    iter = iter + 1;
    feva = feva  + ineva + 2;
end
x = xk;
out = [];
out.value = fk;
out.gradnorm = norm(gk);
out.iter = iter;
out.feva = feva;
end
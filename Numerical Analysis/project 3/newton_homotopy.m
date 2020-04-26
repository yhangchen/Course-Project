function [x,out] = newton_homotopy(f, df, x0, opts)
%{
    Newton's iteration /Broyden algorithm with homotopy

    Parameters
    ----------
    f: 
        objective function
    df:
        jacobian of the objective function
    x0: 
        initial point
    opts.flag:
        '1' for Broyden method, '0' for vanilla method
    opts.num_homo:
        number of homotopy classes, default 10.
    opts.Tol: list
        tolerance
    opts.Lambda: list
        lambdas in the homotopy classes.
    opts.maxIter
        maximum number of iterations, default 1e4

    Returns
    -------
    x: 
        optimal point
    out.iter: 
        number of iterations
    out.err: 
        infinity norm of the residual

%}

% default settings
if ~isfield(opts, 'flag'); opts.flag = 1; end
if ~isfield(opts, 'maxIter'); opts.maxIter = 1e4; end
if ~isfield(opts, 'num_homo'); opts.num_homo = 10; end
if ~isfield(opts, 'Lambda')
    Lambda = zeros(opts.num_homo,1);
    for i = 1:opts.num_homo
        Lambda(i) = i/opts.num_homo;
    end
    opts.Lambda = Lambda; 
end

if ~isfield(opts, 'Tol') 
    Tol = 1e-6 * ones(opts.num_homo,1);
    Tol(end) = eps;
    opts.Tol = Tol; 
end

flag = opts.flag; % whether Broyden
Tol = opts.Tol; % tol of convergence
n = opts.num_homo; % number of homotopy iterations
Lambda = opts.Lambda; % lambda of homotopy method
maxIter = opts.maxIter;
errs = [];
if flag ~= 0 && flag ~= 1
    error('incorrect flag');
end

x = x0; fx0 = f(x0); iter = 0;
for i = 1:n % homotopy iteration
    normf = inf; lambda = Lambda(i); toli = Tol(i);
    h = @(x) f(x) + (lambda-1) * fx0;
    dfx = df(x);
    if flag
        invA = inv(dfx);
    else
        gradf = dfx;
    end
    while normf > toli && iter < maxIter % newton iteration
        iter = iter + 1;
        xp = x;
        if flag
            x = x - invA*h(x);
        else
            x = x - gradf\h(x);
        end
        normf = norm(h(x),inf);
        errs = [errs;norm(f(x),inf)];
        if flag
            g = h(x) - h(xp); y = x - xp;
            denorm = y'*(invA*g);
            if denorm == 0
                continue
            end
            invA = invA - (invA*g - y)*(y'*invA)/denorm;
        else
            gradf = df(x);
        end
    end
end
out = [];
out.iter = iter;
out.err = normf;
out.errs = errs;
end
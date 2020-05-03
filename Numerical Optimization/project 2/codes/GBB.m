function [x, out] = GBB(fun, x0, opts)
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
    opts.alpha0:
        initial stepsize
    opts.sigma1: 
        lower bound of line search
    opts.sigma2: 
        upper bound of line search
    opts.epsilon: 
        bound of the stepsize
    opts.gamma: 
        params in f(x_k-\lambda g_k)\leq \max_{0\leq j\leq \min(k,M)} (f_{k-j})-\gamma\lambda g^\top_k g_k
    opts.M: 
        params f(x_k-\lambda g_k)\leq \max_{0\leq j\leq \min(k,M)} (f_{k-j})-\gamma\lambda g^\top_k g_k
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

if nargin < 3; opts = []; end
% setting default values.
if ~isfield(opts,'linesearch'); opts.linesearch = 'inexact'; end
if ~isfield(opts,'tol'); opts.tol = 1e-6; end
if ~isfield(opts,'alpha0'); opts.alpha0 = 1; end
if ~isfield(opts,'gamma'); opts.gamma = 1e-4; end
if ~isfield(opts,'sigma1'); opts.sigma1 = 0.1; end
if ~isfield(opts,'sigma2'); opts.sigma2 = 0.5; end
if ~isfield(opts,'epsilon'); opts.epsilon = 1e-10; end
if ~isfield(opts,'maxIter'); opts.maxIter = 1e4; end
if ~isfield(opts,'maxinIter'); opts.maxinIter = 10; end
if ~isfield(opts,'M'); opts.M = 10; end
if ~isfield(opts,'debug'); opts.debug = 1; end % print results
x = x0;
n = length(x);
f0 = -inf; 
[f1, g1] = fun(x);
iter = 0; feval = 2;
alpha = opts.alpha0;
fs = []; % encode values of f
xs = {}; values = []; gradnorms = [];
while (norm(g1,inf) > opts.tol*(1+abs(f1))) && abs(f1-f0) > eps
    fs = [fs; f1];
    % step 2
    if alpha <= opts.epsilon || alpha >= 1/opts.epsilon
        % choice of delta
        delta = min(1e5,max(1/norm(g1),1));
        alpha = delta;
    end
    % step 3
    lambda = 1/alpha;
    % step 4/5 nonmonotone line search
    % "Numerical methods for unconstrained optimization and nonlinear equations"
    maxima = max(fs(max(0,end-opts.M)+1:end));
    initer = 0;    % algorithm A6.3.1 in 

    initslope = -lambda*dot(g1, g1);
    lambda_n = 1;
    while true
        x0 = x - lambda_n*lambda*g1;
        [f11,~] = fun(x0);
        if f11 <= maxima + opts.gamma * lambda_n * initslope || initer > opts.maxinIter
            break
        else
            if initer == 0 % first backtrack, quadratic fit
                lambda0 = -initslope / 2 / (f11 - maxima - initslope);
            else % all subsequent backtracks, cubic fit
                ab = 1/(lambda_n - lambda_p)*[1/lambda_n^2, -1/lambda_p^2;...
                    -lambda_p/lambda_n^2, lambda_n/lambda_p^2]*[f11-maxima-lambda_n*initslope;...
                    f11_p - maxima - lambda_p * initslope];
                a = ab(1); b = ab(2);
                disc = b^2 - 3 * a * initslope;
                if abs(a) < eps
                    % cubic is quadratic
                    lambda0 = -initslope/2/b;
                else
                    lambda0 = (disc^0.5-b)/3/a;
                end
            end
            lambda_p = lambda_n;
            f11_p = f11;
            lambda_n = min(max(opts.sigma1*lambda_n, lambda0),opts.sigma2*lambda_n);
        end
        initer = initer + 1;
    end
    lambda = lambda * lambda_n;
    x = x - lambda * g1;
    % step 6
    g0 = g1;
    [f1, g1] = fun(x);
    alpha = -dot(g0,g1-g0)/dot(g0,g0)/lambda;
    feval = feval + 2 + initer;
    iter = iter + 1;
    if opts.debug
        fprintf('iter: '+string(iter) + ', lambda: '+string(lambda) + ', value: ' + string(f1)+', gradnorm: '+ string(norm(g1,inf))+'.\n');
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
out.restart = 0;
out.xs = xs;
out.values = values;
out.gradnorms = gradnorms;
end
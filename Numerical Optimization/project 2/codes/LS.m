function [x, out] = LS(fun, x0, opts)
%{
    LS conjugate gradient methods

    Alg 2/3 in "Efficient Generalized Conjugate Gradient Algorithms"

    Parameters
    ----------
    fun: 
        objective function, with callable method f and g
    x0: 
        initial point
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
if nargin < 3; opts = []; end
if ~isfield(opts,'linesearch'); opts.linesearch = 'inexact'; end
if ~isfield(opts,'tol'); opts.tol = 1e-6; end
if ~isfield(opts,'maxIter'); opts.maxIter = 1e4; end
if ~isfield(opts,'n_restart'); opts.n_restart = opts.maxIter; end % without periodic restart
if ~isfield(opts,'debug'); opts.debug = 1; end

x = x0;
n = length(x);
[f1, g1] = fun(x); alpha = 1;
iter = 0; feval = 2; restart = 0; restart_cond = 1;
nf1 = max(abs(f1),1); ng1 = max(norm(g1),1);
xs = {}; values = []; gradnorms = [];
while 1
    %step 1
    if restart_cond
        d = -g1; initer = 1; restart_cond = 0; restart = restart + 1;
        if opts.debug
            fprintf('Restart.\n');
        end
    end
    % step 2
    alpha0 = alpha;
%     initstep = min(2,-2*f1/dot(g1,alpha*d));
%     Rule.opt(2) = initstep;
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
    iter = iter + 1; initer = initer + 1;
    f0 = f1; g0 = g1;
    [f1, g1] = fun(x);
    feval = feval + info(3) + 2;
    if opts.debug
        fprintf('LS iter: '+string(iter)+', value: ' + string(f1)+', gradnorm: '+ string(norm(g1,inf))+...
            ', directnorm: '+ string(norm(d))+ '.\n');
        xs{iter} = x;
        values = [values;f1];
        gradnorms = [gradnorms;norm(g1,inf)];

    end
    % step 3
    if norm(g1,inf) <= opts.tol*(1+abs(f1)) || iter == opts.maxIter ...
            || abs(f1-f0) <= eps
        break
    end
    % step 4
    if abs(dot(g1,g0)) >= norm(g1)^2/5 ||initer > opts.n_restart% restart every n_restart steps
        restart_cond = 1;
        continue
    end
    % step 5
    % choose delta, gamma, r, eps means the machine epsilon
    r = 1/sqrt(eps);
    delta = sqrt(eps)/norm(d);
    gamma = sqrt(eps)/norm(g0);
    % approximate tk, uk, vk
    [~, g2] = fun(x+delta*d);
    [~, g3] = fun(x+gamma*g1);
    t = dot(d,g2-g1)/delta;
    u = dot(g1,g2-g1)/delta;
    v = dot(g1,g3-g1)/gamma;
    feval = feval + 2;
    % step 6
    if t > 0 && v > 0 && (1-u^2/t/v) >= 1/4/r &&...
            v/dot(g1,g1)/t*dot(d,d) <= r
        % step 7
        d = ((u*dot(g1,d)-t*dot(g1,g1))*g1...
            +(u*dot(g1,g1)-v*dot(g1,d))*d)/(t*v-u^2);
    else
        restart_cond = 1;
        continue
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


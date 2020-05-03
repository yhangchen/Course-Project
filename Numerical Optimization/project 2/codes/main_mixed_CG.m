function [values, gradnorms, iters, fevals, restarts] = main_mixed_CG(choice)
%{
    Test mixed conjugate gradient methods

    method: 
        FR-PRP / HS

    Parameters
    ----------
    choice = 1,2,3,4
        Specify test functions.

    Returns
    -------
    values, gradnorms, iters, fevals, restarts:
        matrices encoding the results.
%}

addpath('linesearch')
addpath('test')

ns = [100,1000,10000];
if choice == 4
    ns = [100,1024,10000];
end
values = zeros(2,length(ns));
gradnorms = zeros(2,length(ns));
iters = zeros(2,length(ns));
fevals = zeros(2,length(ns));
restarts = zeros(2,length(ns));
times = zeros(2,length(ns));


for i = 1 : length(ns)
    n = ns(i);
    % set test functions and default init
    if choice == 1
        fun = @(x) test1(x,n);
        x0 = ones(n, 1) / n;
    elseif choice == 2
        fun = @(x) test2(x,n);
        x0 = repmat([3,-1,0,3],1,n/4)';
    elseif choice == 3
        fun = @(x) test3(x,n);
        x0 = ones(n, 1)/n^2;
    elseif choice == 4
        n0 = sqrt(n);
        if abs(n0^2 - n) > eps
            error('n should be a square of an integer.')
        end
        fun = @(x) test4(x,floor(n0));
        x0 = (sin( [1 : n0^2].^2 )' / 5);
    else
        error('incorrect test function');
    end
    tic
    [x1, out1] = CG(fun, x0, 'FR-PRP');
    t1 = toc;
    tic
    [x2, out2] = LS(fun, x0);
    t2 = toc;
    values(1,i) = out1.value; values(2,i) = out2.value; 
    gradnorms(1,i) = out1.gradnorm; gradnorms(2,i) = out2.gradnorm; 
    iters(1,i) = out1.iter; iters(2,i) = out2.iter;
    fevals(1,i) = out1.feval; fevals(2,i) = out2.feval; 
    restarts(1,i) = out1.restart; restarts(2,i) = out2.restart;
    times(1,i) = t1; times(2,i) = t2;
end
values = values';
gradnorms = gradnorms';
iters = iters';
fevals = fevals';
restarts = restarts'-1;
times = times';
end
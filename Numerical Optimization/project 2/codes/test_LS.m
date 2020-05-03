function results = test_LS(choice)
%{
    Test non-linear conjugate gradient methods

    method: 
        'FR' for FR, 'PRP' for PRP, 'PRP+' for PRP+, 'HS' for HS,
             'CD' for conjugate descent, 'DY' for Dai-Yuan.
         FR-PRP / HS
        
         Global BB methods

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

for i = 1 : length(ns)
    n = ns(i);
    if choice == 1
        fun = @(x) test1(x,n);
        x0 = ones(n, 1) / n;
    elseif choice == 2
        fun = @(x) test2(x,n);
        x0 = repmat([3,-1,0,3],1,n/4)';
    elseif choice == 3
        fun = @(x) test3(x,n);
        x0 = ones(n, 1);
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
    [x1, out1] = LS(fun, x0);
    t1 = toc;
    values(1,i) = out1.value;
    gradnorms(1,i) = out1.gradnorm; 
    iters(1,i) = out1.iter; 
    fevals(1,i) = out1.feval; 
    restarts(1,i) = out1.restart; 
    times(1,i) = t1;
end
values = values';
gradnorms = gradnorms';
iters = iters';
fevals = fevals';
restarts = restarts'-1;
times = times';

results = [values;gradnorms;iters;fevals;restarts;times];


end
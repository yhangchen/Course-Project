function results = main_test(choice)
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
values = zeros(9,length(ns));
gradnorms = zeros(9,length(ns));
iters = zeros(9,length(ns));
fevals = zeros(9,length(ns));
restarts = zeros(9,length(ns));
times = zeros(9,length(ns));

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
    [x1, out1] = CG(fun, x0, 'FR');
    t1 = toc;
    tic
    [x2, out2] = CG(fun, x0, 'PRP');
    t2 = toc; 
    tic
    [x3, out3] = CG(fun, x0, 'PRP+');
    t3 = toc;
    tic
    [x4, out4] = CG(fun, x0, 'CD');
    t4 = toc;
    tic
    [x5, out5] = CG(fun, x0, 'DY');
    t5 = toc;
    tic
    [x6, out6] = CG(fun, x0, 'HS');
    t6 = toc;
    tic
    [x7, out7] = CG(fun, x0, 'FR-PRP');
    t7 = toc;
    tic
    [x8, out8] = LS(fun, x0);
    t8 = toc;
    tic
    [x9, out9] = GBB(fun, x0);
    t9 = toc;
    values(1,i) = out1.value; values(2,i) = out2.value; values(3,i) = out3.value; values(4,i) = out4.value; values(5,i) = out5.value; values(6,i) = out6.value;values(7,i) = out7.value; values(8,i) = out8.value;values(9,i) = out9.value;
    gradnorms(1,i) = out1.gradnorm; gradnorms(2,i) = out2.gradnorm; gradnorms(3,i) = out3.gradnorm; gradnorms(4,i) = out4.gradnorm; gradnorms(5,i) = out5.gradnorm; gradnorms(6,i) = out6.gradnorm;gradnorms(7,i) = out7.gradnorm; gradnorms(8,i) = out8.gradnorm;gradnorms(9,i) = out9.gradnorm;
    iters(1,i) = out1.iter; iters(2,i) = out2.iter; iters(3,i) = out3.iter; iters(4,i) = out4.iter; iters(5,i) = out5.iter; iters(6,i) = out6.iter;iters(7,i) = out7.iter; iters(8,i) = out8.iter; iters(9,i) = out9.iter;
    fevals(1,i) = out1.feval; fevals(2,i) = out2.feval; fevals(3,i) = out3.feval; fevals(4,i) = out4.feval; fevals(5,i) = out5.feval; fevals(6,i) = out6.feval;fevals(7,i) = out7.feval; fevals(8,i) = out8.feval; fevals(9,i) = out9.feval;
    restarts(1,i) = out1.restart; restarts(2,i) = out2.restart; restarts(3,i) = out3.restart; restarts(4,i) = out4.restart; restarts(5,i) = out5.restart; restarts(6,i) = out6.restart;restarts(7,i) = out7.restart; restarts(8,i) = out8.restart; restarts(9,i) = out9.restart;
    times(1,i) = t1; times(2,i) = t2;times(3,i) = t3;times(4,i) = t4;times(5,i) = t5;times(6,i) = t6; times(7,i) = t7;times(8,i) = t8;times(9,i) = t9;
end


values = values';
gradnorms = gradnorms';
iters = iters';
fevals = fevals';
restarts = restarts'-1;
times = times';

results = [values;gradnorms;iters;fevals;restarts;times];

end
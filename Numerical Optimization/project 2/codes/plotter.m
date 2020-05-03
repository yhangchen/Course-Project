function plotter(choice)
%{
    Test non-linear conjugate gradient methods

    method: 
        'FR' for FR, 'PRP' for PRP, 'PRP+' for PRP+, 'HS' for HS,
             'CD' for conjugate descent, 'DY' for Dai-Yuan.

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
its = [250,5000,10000]; %[200,200,200]
values = zeros(6,length(ns));
gradnorms = zeros(6,length(ns));
iters = zeros(6,length(ns));
fevals = zeros(6,length(ns));
restarts = zeros(6,length(ns));

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
    opts.maxIter = its(i);
    [x1, out1] = CG(fun, x0, 'FR',opts);
        
    [x2, out2] = CG(fun, x0, 'PRP',opts);
    
    [x3, out3] = CG(fun, x0, 'PRP+',opts);
    
    [x4, out4] = CG(fun, x0, 'CD',opts);
    
    [x5, out5] = CG(fun, x0, 'DY',opts);
    
    [x6, out6] = CG(fun, x0, 'HS',opts);
    
    [x7, out7] = CG(fun, x0, 'FR-PRP',opts);
    
    [x8, out8] = LS(fun, x0,opts);
    
    [x8, out9] = GBB(fun, x0,opts);

    fig = figure;
    name = 'report\'+string(choice)+'_'+string(n)+'_value.eps';
    semilogy(out1.values,'Linewidth',2); hold on
    semilogy(out2.values,'Linewidth',2);
    semilogy(out3.values,'Linewidth',2);
    semilogy(out4.values,'Linewidth',2);
    semilogy(out5.values,'Linewidth',2);
    semilogy(out6.values,'Linewidth',2);
    semilogy(out7.values,'Linewidth',2);
    semilogy(out8.values,'Linewidth',2);
    semilogy(out9.values,'Linewidth',2);
    legend('FR','PRP','PRP+','CD','DY','HS','FR-PRP','LS','GBB','interpreter','latex','fontsize',10);
    xlabel('iteration','interpreter','latex','fontsize',12);
    ylabel('$f_k$','interpreter','latex','fontsize',12);
    xlim([1,its(i)]);
    grid on
    print(fig,name,'-depsc2')
    hold off

    fig = figure;
    name = 'report\'+string(choice)+'_'+string(n)+'_grad.eps';
    semilogy(out1.gradnorms,'Linewidth',2); hold on
    semilogy(out2.gradnorms,'Linewidth',2);
    semilogy(out3.gradnorms,'Linewidth',2);
    semilogy(out4.gradnorms,'Linewidth',2);
    semilogy(out5.gradnorms,'Linewidth',2);
    semilogy(out6.gradnorms,'Linewidth',2);
    semilogy(out7.gradnorms,'Linewidth',2);
    semilogy(out8.gradnorms,'Linewidth',2);
    semilogy(out9.gradnorms,'Linewidth',2);
    legend('FR','PRP','PRP+','CD','DY','HS','FR-PRP','LS','BB','interpreter','latex','fontsize',10);
    xlabel('iteration','interpreter','latex','fontsize',12);
    ylabel('$\|g_k\|_\infty$','interpreter','latex','fontsize',12);
    xlim([1,its(i)]);
    grid on
    print(fig,name,'-depsc2')
    hold off
end


values = values';
gradnorms = gradnorms';
iters = iters';
fevals = fevals';
restarts = restarts'-1;
end
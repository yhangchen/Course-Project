% min 0.5 ||Ax-b||_2^2 + mu*||x||_1,2
clear all
seed = 97006855;
ss = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(ss);
n = 512;
m = 256;
k = round(n*0.1); l = 2;
A = randn(m,n);
p = randperm(n); p = p(1:k);
u = zeros(n,l); u(p,:) = randn(k,l);
b = A*u;
mu = 1e-2;


x0 = rand(n,l);


errfun = @(x1, x2) norm(x1-x2,'fro')/(1+norm(x1,'fro'));
errfun_exact = @(x1) norm(x1-u,'fro')/(1+norm(u,'fro'));
sparsity = @(x) sum(abs(x(:))>1e-6*max(abs(x(:))))/(n*l);

% cvx calling mosek
opts1 = [];
tic;
[x1, iter1, out1] = gl_cvx_mosek(x0, A, b, mu, opts1);
t1 = toc;

% cvx calling gurobi
opts2 = [];
tic;
[x2, iter2, out2] = gl_cvx_gurobi(x0, A, b, mu, opts2);
t2 = toc;

% call mosek directly
opts3 = [];
tic;
[x3, iter3, out3] = gl_mosek(x0, A, b, mu, opts3);
t3 = toc;

% call gurobi directly
opts4 = [];
tic;
[x4, iter4, out4] = gl_gurobi(x0, A, b, mu, opts4);
t4 = toc;

opts5 = [];
tic;
[x5, iter5, out5] = gl_SGD_primal(x0, A, b, mu, opts5);
t5 = toc;

opts6 = [];
tic;
[x6, iter6, out6] = gl_GD_primal(x0, A, b, mu, opts6);
t6 = toc;

opts7 = [];
tic;
[x7, iter7, out7] = gl_FGD_primal(x0, A, b, mu, opts7);
t7 = toc;

opts8 = [];
tic;
[x8, iter8, out8] = gl_ProxGD_primal(x0, A, b, mu, opts8);
t8 = toc;

opts9 = [];
tic;
[x9, iter9, out9] = gl_FProxGD_primal(x0, A, b, mu, opts9);
t9 = toc;

opts10 = [];
tic;
[x10, iter10, out10] = gl_ALM_dual(x0, A, b, mu, opts10);
t10 = toc;


opts11 = [];
tic;
[x11, iter11, out11] = gl_ADMM_dual(x0, A, b, mu, opts11);
t11 = toc;

opts12 = [];
tic;
[x12, iter12, out12] = gl_ADMM_primal(x0, A, b, mu, opts12);
t12 = toc;


fig = figure;
hold off
subplot(13,1,1)
plot(u)
hold on
legend('solution')
subplot(13,1,2)
plot(x1)
legend('cvx-mosek')
subplot(13,1,3)
plot(x2)
legend('cvx-gurobi')
subplot(13,1,4)
plot(x3)
legend('mosek')
subplot(13,1,5)
plot(x4)
legend('gurobi');
subplot(13,1,6)
plot(x5)
legend('subgrad')
subplot(13,1,7)
plot(x6)
legend('grad smooth');
subplot(13,1,8)
plot(x7)
legend('fast smooth');
subplot(13,1,9)
plot(x8)
legend('proxgrad');
subplot(13,1,10)
plot(x9)
legend('fast proxgrad');
subplot(13,1,11)
plot(x10)
legend('alm dual');
subplot(13,1,12)
plot(x11)
legend('admm dual');
subplot(13,1,13)
plot(x12)
legend('admm lprimal');


hold off

results = [
    t1, out1.val, errfun(x1, x1), errfun_exact(x1), sparsity(x1), iter1
    t2, out2.val, errfun(x1, x2), errfun_exact(x2), sparsity(x2), iter2
    t3, out3.val, errfun(x1, x3), errfun_exact(x3), sparsity(x3), iter3
    t4, out4.val, errfun(x1, x4), errfun_exact(x4), sparsity(x4), iter4
    t5, out5.val, errfun(x1, x5), errfun_exact(x5), sparsity(x5), iter5
    t6, out6.val, errfun(x1, x6), errfun_exact(x6), sparsity(x6), iter6
    t7, out7.val, errfun(x1, x7), errfun_exact(x7), sparsity(x7), iter7
    t8, out8.val, errfun(x1, x8), errfun_exact(x8), sparsity(x8), iter8
    t9, out9.val, errfun(x1, x9), errfun_exact(x9), sparsity(x9), iter9
    t10, out10.val, errfun(x1, x10), errfun_exact(x10), sparsity(x10), iter10
    t11, out11.val, errfun(x1, x11), errfun_exact(x11), sparsity(x11), iter11
    t12, out12.val, errfun(x1, x12), errfun_exact(x12), sparsity(x12), iter12
    ]


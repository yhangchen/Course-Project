% function Test_group_lasso

% min 0.5 * ||A * x - b||_2^2 + mu * ||x||_{1,2}

% generate data
m=1024;n=512;
A = load("../results/Frob_L_0_BB_A.csv");
b = load("../results/Frob_L_0_BB_b.csv");
u = load("../results/Frob_L_0_BB_exact_x.csv");
mu = 1e-2;
l = 1; L = 4;
x0 = load("../results/Frob_L_0_BB_x0.csv");
opt_lst = {"l12", "l21",{"elastic",0.5},"linf"};


% change to the following for matrix
% A = load("../results/Frob_L_12_BB_A.csv");
% b = load("../results/Frob_L_12_BB_b.csv");
% u = load("../results/Frob_L_12_BB_exact_x.csv");
% mu = 1e-2;
% l = 4; L = 3;
% x0 = load("../results/Frob_L_12_BB_x0.csv");
% opt_lst = {"l12", "l21",{"elastic",0.5}};


errfun = @(x1, x2) norm(x1 - x2, 'fro') / (1 + norm(x1,'fro'));
errfun_exact = @(x, l) norm(x(:, 1: l) - u(:, 1: l), 'fro') /...
    (1 + norm(u(:, 1: l),'fro'));
sparisity = @(x, l) sum(abs(x(:)) > 1E-6 * max(abs(x(:)))) /(n*l);
%, "linf", "TV1D", "TV2D", "nuclear", {"ind_nuclear", 1}
x_lst = cell(1, L);
iter_lst = zeros(1, L);
out_lst = zeros(1, L);
penal_lst = zeros(1, L);

t_lst = zeros(1, L);
% l_lst = [l, l, 1, 1, l, l, l];
l_lst = [l,l,l,l];
% cvx calling mosek
for i = 1: L
    tic;
    [x_lst{i}, iter_lst(i), out_lst(i), penal_lst(i)] =...
        gl_cvx_mosek(x0(:, 1: l_lst(i)), A, b(:, 1: l_lst(i)), mu, opt_lst{i});
    t_lst(i) = toc;
end

% print comparison results with cvx-call-mosek
for i = 1: L
    opt = opt_lst{i};
    mode = opt{1};
    temp_x = x_lst{i};
    fprintf("CVX-Mosek-" + mode +...
        ": cpu: %5.2f, iter: %5d, optval: %6.5E, penalty val: %6.5E, sparisity: %4.3f,err-to-exact: %3.2E.\n",...
        t_lst(i), iter_lst(i), out_lst(i), penal_lst(i), sparisity(...
        temp_x(:, 1: l_lst(i)), l_lst(i)),...
        errfun_exact(temp_x(:, 1: l_lst(i)), l_lst(i)));
end
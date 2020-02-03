function [output] = inexact_Uzawa_builtin(choice)
% implementation of inexact Uzawa based on V-cycle.
    n = choice.n;
	L = choice.L;
    if ~(n>=2) || n <= L
        error('Wrong size')
    end
    tic
    alpha = choice.alpha;
    N = 2^n;
    max_iter_out = choice.max_iter_out;
    iter = 0;
    eps = choice.eps;
	U = zeros(N*(N-1),1); V = zeros(N*(N-1),1); P = zeros(N^2,1);
    F = mat_F(N); G = mat_G(N); A = N^2*mat_A(N); B1 = N * mat_B1(N); B2 = N * mat_B2(N);
    F_r = F; G_r = G; rF = F; rG = G;
    init = norm([F;G]);
    while norm([F_r;G_r]) > init * eps   && iter < max_iter_out
        U = A\rF; V = A\rG;
        P = P + alpha * (B1'*U + B2'*V);
        rF = F - B1 * P; F_r = rF - A * U;
        rG = G - B2 * P; G_r = rG - A * V;
        iter = iter + 1;
    end
 
	time = toc;
	output = [];
	output.time = time;
	output.iter = iter;
	output.U = U;
	output.V = V;
	output.P = P;
	output.res = norm([F_r;G_r])/init;
    output.u_err = norm([U;V]-mat_U(N))/N;
    output.p_err = norm(P-mat_P(N))/N;
end



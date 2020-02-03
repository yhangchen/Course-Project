function [output] = inexact_Uzawa_CG(choice)
% implementation of inexact Uzawa based on PCG.   
    n = choice.n;
    tau = choice.tau;
    tic
    alpha = choice.alpha;
    N = 2^n;
    max_iter_in = choice.max_iter_in;
    max_iter_out = choice.max_iter_out;
    iter = 0;
	U = zeros(N*(N-1),1); V = zeros(N*(N-1),1); P = zeros(N^2,1);
    F = mat_F(N); G = mat_G(N); A = N^2*mat_A(N); B1 = N * mat_B1(N); B2 = N * mat_B2(N);
    F_r = F; G_r = G;
    eps = choice.eps;
    D = spdiags(diag(A),0,N*(N-1),N*(N-1));
    init = norm([F_r;G_r]);
    while norm([F_r;G_r]) > init * eps && iter < max_iter_out
		[U,~,~,~,~] = pcg(A, F - B1 * P, tau, max_iter_in, D); [V,~,~,~,~] = pcg(A, G - B2 * P, tau, max_iter_in, D);
% U = pcg(A, F - B1 * P, tau, max_iter_in, D); V = pcg(A, G - B2 * P, tau, max_iter_in, D);
        P = P + alpha * (B1'*U + B2'*V);
        F_r = F - B1 * P - A * U;
        G_r = G - B2 * P - A * V;
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

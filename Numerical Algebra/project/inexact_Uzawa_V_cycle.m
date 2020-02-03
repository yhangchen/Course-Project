function [output] = inexact_Uzawa_V_cycle(choice)
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
    eps = choice.eps;
    A_store = {}; Rt_u = {};
    for k=1:L+1
		Nk = 2^(n-k+1);
		A_store{k} = Nk^2 * mat_A(Nk);
        if k<=L
			Rt_u{k} = restrict_u(Nk);
        end
    end
    
    choice.A_store = A_store;
    choice.Rt_u = Rt_u;
    
    iter = 0;
	U = zeros(N*(N-1),1); V = zeros(N*(N-1),1); P = zeros(N^2,1);
    F = mat_F(N); G = mat_G(N); B1 = N * mat_B1(N); B2 = N * mat_B2(N);
    F_r = F; G_r = G; rF = F; rG = G;
    init = norm([F;G]);
    while norm([F_r;G_r]) > init * eps   && iter < max_iter_out
        choice.F = rF; 
        choice.G = rG; 
        choice.U = U;
        choice.V = V;
        choice.tau0 = choice.tau;
        [U, V, iter_in] = V_cycle_Uzawa_inner(choice);
		P = P + alpha * (B1'*U + B2'*V);
        rF = F - B1 * P; F_r = rF - A_store{1} * U;
        rG = G - B2 * P; G_r = rG - A_store{1} * V;
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


function [Uk,Vk,iter] = V_cycle_Uzawa_inner(choice)

	v1 = choice.v1;
	v2 = choice.v2;
    tau0 = choice.tau0;
    U0 = choice.U;
    V0 = choice.V;
    F0 = choice.F;
    G0 = choice.G;
    A_store = choice.A_store;
    Rt_u = choice.Rt_u;
    n = choice.n; N = 2^n; L = choice.L;
    max_iter_in = choice.max_iter_in;
	U_store = {}; V_store = {}; F_store = {}; G_store = {};
    
    Uk = Gauss_Seidel(A_store{1},F0,U0,v1);
    Vk = Gauss_Seidel(A_store{1},G0,V0,v1);

	U_store{1}=Uk;V_store{1}=Vk;F_store{1}=F0;G_store{1}=G0;
	F_r = F0-A_store{1}*Uk;G_r = G0-A_store{1}*Vk;

	iter = 0; 
    while norm([F_r;G_r]) > N^2*tau0 && iter < max_iter_in
		Fk = Rt_u{1}*F_r; Gk = Rt_u{1}*G_r;
		k = 1;
        while k<L
			Nk = 2^(n-k);
			Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1);
            Uk = Gauss_Seidel(A_store{k+1},Fk,Uk,v1);
            Vk = Gauss_Seidel(A_store{k+1},Gk,Vk,v1);
			U_store{k+1}=Uk;V_store{k+1}=Vk;
			F_store{k+1}=Fk;G_store{k+1}=Gk;
			Fk = Rt_u{k+1}*(Fk-A_store{k+1}*Uk);
			Gk = Rt_u{k+1}*(Gk-A_store{k+1}*Vk);
			k = k+1;
        end        
		Nk = 2^(n-k);        
		Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1);
        Uk = Gauss_Seidel(A_store{k+1},Fk,Uk,v1);
        Vk = Gauss_Seidel(A_store{k+1},Gk,Vk,v1);     
        while k>0
			Uk = U_store{k}+4*(Rt_u{k})'*Uk;
			Vk = V_store{k}+4*(Rt_u{k})'*Vk;
            Uk = Gauss_Seidel(A_store{k},F_store{k},Uk,v2);
            Vk = Gauss_Seidel(A_store{k},G_store{k},Vk,v2);
			k = k-1;
        end  
		F_r = F0-A_store{1}*Uk; G_r = G0-A_store{1}*Vk;
		U_store{1} = Uk; V_store{1} = Vk;
		iter = iter+1;
    end
end

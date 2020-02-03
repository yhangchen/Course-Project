function [output] = V_cycle(choice)
    n = choice.n;
	L = choice.L;
    
    if ~(n>=2) || n <= L
        error('Wrong size')
    end
    
	smoother = choice.smoother;
	v1 = choice.v1;
	v2 = choice.v2;
	eps = choice.eps;
    if strcmp(smoother,'exact')||strcmp(smoother,'inexact')
        alpha = choice.alpha;
    end
    N = 2^n;
    max_iter = choice.max_iter;
    % construct the matrix in the V-cycle process.
    A_store = {}; B1_store = {}; B2_store = {};
    Rt_u = {}; Rt_p = {};
    U_store = {}; V_store = {}; P_store = {}; F_store = {}; G_store = {}; D_store = {};
    if strcmp(smoother, '1')
        BtB_store = {};
    end
    
    tic
    for k=1:L+1
		Nk = 2^(n-k+1);
		A_store{k} = Nk^2 * mat_A(Nk); B1_store{k} = Nk * mat_B1(Nk); B2_store{k} = Nk * mat_B2(Nk);
        if strcmp(smoother, '1')
			BtB_store{k} = Nk^2 * mat_BtB(Nk);
        end
        if k<=L
			Rt_u{k} = restrict_u(Nk);
			Rt_p{k} = restrict_p(Nk);
        end
    end
    
    F0 = mat_F(N); G0 = mat_G(N);
	Uk = zeros(N*(N-1),1); Vk = zeros(N*(N-1),1); Pk = zeros(N^2,1); Dk = zeros(N^2,1);

    if strcmp(smoother,'1')
		[Uk,Vk,Pk] = DGS_1(Uk,Vk,Pk,F0,G0,A_store{1},B1_store{1},B2_store{1},BtB_store{1}, Dk, N, v1);
    elseif strcmp(smoother,'2')
        [Uk,Vk,Pk] = DGS_2(Uk,Vk,Pk,F0,G0,A_store{1},B1_store{1},B2_store{1}, Dk, N, v1);
    elseif strcmp(smoother,'exact')
		[Uk,Vk,Pk] = Uzawa(Uk,Vk,Pk,F0,G0,A_store{1},B1_store{1},B2_store{1}, Dk, N, v1, alpha);
    elseif strcmp(smoother,'inexact')
        [Uk,Vk,Pk] = inexact_Uzawa(Uk,Vk,Pk,F0,G0,A_store{1},B1_store{1},B2_store{1}, Dk, N, v1, alpha);
    end 

    iter = 0;
    
	U_store{1}=Uk;V_store{1}=Vk;P_store{1}=Pk;F_store{1}=F0;G_store{1}=G0;D_store{1}=Dk;
	F_r = F0-A_store{1}*Uk-B1_store{1}*Pk; G_r = G0-A_store{1}*Vk-B2_store{1}*Pk; D_r = -B1_store{1}'*Uk-B2_store{1}'*Vk;
    init = norm([F0;G0]);
    while norm([F_r;G_r;D_r],'fro') > init*eps && iter < max_iter
        k = 1;
		Fk = Rt_u{1}*F_r; Gk = Rt_u{1}*G_r; Dk = Rt_p{1}*D_r;
        while k<L
			Nk = 2^(n-k);
			Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1); Pk = zeros(Nk^2,1);
            if strcmp(smoother,'1')
                [Uk,Vk,Pk] = DGS_1(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1},BtB_store{k+1}, Dk, Nk, v1);
            elseif strcmp(smoother,'2')
                [Uk,Vk,Pk] = DGS_2(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1}, Dk, Nk, v1);
            elseif strcmp(smoother,'exact')
                [Uk,Vk,Pk] = Uzawa(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1}, Dk, Nk, v1, alpha);
            elseif strcmp(smoother,'inexact')
                [Uk,Vk,Pk] = inexact_Uzawa(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1}, Dk, Nk, v1, alpha);
            end 
			U_store{k+1}=Uk;V_store{k+1}=Vk;P_store{k+1}=Pk;F_store{k+1}=Fk;G_store{k+1}=Gk;D_store{k+1}=Dk;
			Fk = Rt_u{k+1}*(Fk-A_store{k+1}*Uk-B1_store{k+1}*Pk);
			Gk = Rt_u{k+1}*(Gk-A_store{k+1}*Vk-B2_store{k+1}*Pk);
			Dk = Rt_p{k+1}*(Dk-B1_store{k+1}'*Uk-B2_store{k+1}'*Vk);
			k = k+1;
        end
		Nk = 2^(n-k);
		Uk = zeros(Nk*(Nk-1),1); Vk = zeros(Nk*(Nk-1),1); Pk = zeros(Nk^2,1);
        if strcmp(smoother,'1')
            [Uk,Vk,Pk] = DGS_1(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1},BtB_store{k+1}, Dk, Nk, v1);
        elseif strcmp(smoother,'2')
            [Uk,Vk,Pk] = DGS_2(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1}, Dk, Nk, v1);
        elseif strcmp(smoother,'exact')
            [Uk,Vk,Pk] = Uzawa(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1},Dk, Nk, v1, alpha);
        elseif strcmp(smoother,'inexact')
            [Uk,Vk,Pk] = inexact_Uzawa(Uk,Vk,Pk,Fk,Gk,A_store{k+1},B1_store{k+1},B2_store{k+1},Dk, Nk, v1, alpha);
        end 
        while k>0
            Nk = 2^(n-k+1);
			Uk = U_store{k}+4*(Rt_u{k})'*Uk;
			Vk = V_store{k}+4*(Rt_u{k})'*Vk;
			Pk = P_store{k}+4*(Rt_p{k})'*Pk;
			
            if strcmp(smoother,'1')
                [Uk,Vk,Pk] = DGS_1(Uk,Vk,Pk,F_store{k},G_store{k},A_store{k},B1_store{k},B2_store{k},BtB_store{k}, D_store{k}, Nk, v2);
            elseif strcmp(smoother,'2')
                [Uk,Vk,Pk] = DGS_2(Uk,Vk,Pk,F_store{k},G_store{k},A_store{k},B1_store{k},B2_store{k}, D_store{k}, Nk, v2);
            elseif strcmp(smoother,'exact')
                [Uk,Vk,Pk] = Uzawa(Uk,Vk,Pk,F_store{k},G_store{k},A_store{k},B1_store{k},B2_store{k}, D_store{k}, Nk, v2, alpha);
            elseif strcmp(smoother,'inexact')
                [Uk,Vk,Pk] = inexact_Uzawa(Uk,Vk,Pk,F_store{k},G_store{k},A_store{k},B1_store{k},B2_store{k}, D_store{k}, Nk, v2, alpha);
            end             
			k = k-1;
        end
		U_store{1} = Uk;V_store{1} = Vk;P_store{1} = Pk;
		F_r = F0-A_store{1}*Uk-B1_store{1}*Pk;
		G_r = G0-A_store{1}*Vk-B2_store{1}*Pk;
		D_r = -B1_store{1}'*Uk-B2_store{1}'*Vk;
		iter = iter+1;
    end
    
	time = toc;
	output = [];
	output.time = time;
	output.iter = iter;
	output.U = Uk;
	output.V = Vk;
	output.P = Pk;
	output.res = norm([F_r;G_r;D_r])/init;
    output.u_err = norm([Uk;Vk]-mat_U(N))/N;
    output.p_err = norm(Pk-mat_P(N))/N;
end

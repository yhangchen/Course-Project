function Inexact_Uzawa_CG_test
N = [6,7,8];
Alpha = [0.5,1,1.5,2];
Tau = [16*1e-9,4*1e-9,1e-9];


for i = 1:length(N)
for k1 = 1:length(Alpha)
for k2 = 1:length(Tau)
choice=[];
choice.alpha = Alpha(k1);
choice.max_iter_out = 100;
choice.max_iter_in = 1000;
choice.eps = 1e-8;
choice.tau = Tau(k2);
choice.n = N(i);
out = inexact_Uzawa_CG(choice);
fprintf('%d, %.1f, %.2e, %d, %.2f,%.2e, %.2e, %.2e\n',2^N(i),Alpha(k1),Tau(k2),out.iter,out.time,out.res,out.u_err,out.p_err)
end
end
end
end


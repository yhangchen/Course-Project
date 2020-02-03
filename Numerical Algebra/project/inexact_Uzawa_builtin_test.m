function inexact_Uzawa_builtin_test
N = [11];
NmL = [2];
Alpha = [1];

for i = 1:length(N)
for k = 1:length(NmL)
for k1 = 1:length(Alpha)
if NmL(k)<N(i)
choice=[];
choice.alpha = Alpha(k1);
choice.max_iter_out = 100;
choice.max_iter_in = 100;
choice.L = N(i)-NmL(k);
choice.eps = 1e-8;
choice.n = N(i);
out = inexact_Uzawa_builtin(choice);
fprintf('%d, %d, %d, %d, %.2f, %.2e, %.2e, %.2e\n',2^N(i),Alpha(k1),NmL(k),out.iter,out.time,out.res,out.u_err,out.p_err)
end
end
end
end
end
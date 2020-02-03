function inexact_Uzawa_V_cycle_test
% N = [6,7,8,9];
% V = [2,4,6];
% NmL = [1,2];
% Alpha = [1]; 
% Tau = [16*1e-8,4*1e-8,1e-8];

% N = [10,11];
% V = [1,2,3,4];
% NmL = [1];
% Alpha = [1];
% Tau = [1e-8,0.5*1e-8,1e-9];


for i = 1:length(N)
for k = 1:length(NmL)
for j = 1:length(V)
for k1 = 1:length(Alpha)
for k2 = 1:length(Tau)
if NmL(k)<N(i)
choice=[];
choice.alpha = Alpha(k1);
choice.max_iter_out = 10;
choice.max_iter_in = 100;
choice.L = N(i)-NmL(k);
choice.v1 = V(j);
choice.v2 = V(j);
choice.eps = 1e-8;
choice.tau = Tau(k2);
choice.n = N(i);
out = inexact_Uzawa_V_cycle(choice);
fprintf('%d, %d, %d, %d, %.1e, %d, %.2f, %.2e, %.2e, %.2e\n',2^N(i),Alpha(k1), NmL(k),V(j),Tau(k2),out.iter,out.time,out.res,out.u_err,out.p_err)
end
end
end
end
end
end
end
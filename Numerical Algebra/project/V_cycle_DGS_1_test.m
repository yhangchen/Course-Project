function V_cycle_DGS_1_test
% N = [6,7,8,9];
% V = [2,4,8,16,32];
% NmL = [1,2,3];
% N = [10,11];
% V = [2,4,8];
% NmL = [1,2];
% 


for i = 1:length(N)
for k = 1:length(NmL)
for j = 1:length(V)
if NmL(k)<N(i)
choice=[];
choice.smoother = '1';
choice.max_iter = 100;
choice.L = N(i)-NmL(k);
choice.v1 = V(j);
choice.v2 = V(j);
choice.eps = 1e-8;
choice.n = N(i);
out = V_cycle(choice);
fprintf('%d, %d, %d, %d, %.2f,%.2e, %.2e,%.2e\n',2^N(i),NmL(k),V(j),out.iter,out.time,out.res,out.u_err,out.p_err)
end
end
end
end
end